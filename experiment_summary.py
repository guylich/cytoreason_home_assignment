import requests
import pandas as pd
import xmltodict
import json

ENTREZ_SEARCH_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db={db_name}&term={search_term}&retmode=json&retmax=1000'
ENTREZ_SUMMARY_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db={db_name}&id={ids}&retmode=json'
MICROARRAY_PATH = 'results/microarray/{gse_id}_experiment_summary.csv'
RNASEQ_PATH = 'results/rnaseq/{gse_id}_rnaseq.csv'

GSE_DB_NAME = 'gds'
SRA_DB_NAME = 'sra'

def entrez_search(db_name: str, term: str) -> list:
    """
    Queries nbci DB using the Entrez search method.
    Returns the IDs found.
    """
    response = requests.get(ENTREZ_SEARCH_URL.format(db_name=db_name, search_term=term))
    if response.ok:
        return response.json()['esearchresult']['idlist']
    return []


def entrez_summary(db_name: str, ids: list) -> dict:
    """
    Returns a summary of the IDs using the Entrez summary method.
    """
    response=requests.get(ENTREZ_SUMMARY_URL.format(db_name=db_name, ids=','.join(ids)))
    if response.ok:
        return response.json()['result']
    return {}


def sum_GSE(gse_id):
    """
    Recieves a GSE ID and returns its summary
    """
    query ='{gse_id}+AND+gse[ETYP]'.format(gse_id=gse_id)  # returns only the GSE (no GPL or samples)
    ids = entrez_search(GSE_DB_NAME, query)
    summary = entrez_summary(GSE_DB_NAME, ids)
    return summary


def parse_gse_summary(gse_summray: dict) -> pd.DataFrame:
    """
    Extract the relevant fields from a GSE experiment summary into a dataframe
    """
    parsed = list()
    for uid in gse_summray['uids']: 
        parsed.append(dict(
            uid=uid,
            gpl=gse_summray[uid]['gpl'],
            suppfile=gse_summray[uid]['suppfile'],
            ftplink=gse_summray[uid]['ftplink'],
        ))
    return pd.DataFrame(parsed)


def find_relevant_srp(gse_summary: dict) -> list:
    """
    Receives a summary of a GSE experiment and returns SRP IDs that are related to it.
    This is a way to determine if a GSE experiment has any RNAseq data that has to be processed.
    """
    extrelations = []
    for uid in gse_summary['uids']:
        extrelations.extend(gse_summary[uid]['extrelations'])
    srp_ids = [relation['targetobject'] for relation in extrelations if relation['relationtype'] == 'SRA']
    return srp_ids


def parse_sra_summary(sra_summary: dict) -> pd.DataFrame:
    """
    Extracts relevant fields from an SRA summary into a dataframe.
    The field extraction is quite ugly at the moment. With more time and understanding of the data it can be improved.
    """
    parsed = []
    for uid in sra_summary['uids']:
        sra_summary[uid]['expxml'] = xmltodict.parse(f'<root>{sra_summary[uid]["expxml"]}</root>')['root'] # I need to wrap this XML bit with a tag for the parser to work
        sra_summary[uid]['runs'] = xmltodict.parse(sra_summary[uid]["runs"])
        parsed.append(dict(
            uid=uid,
            run_id=sra_summary[uid]['runs']['Run']['@acc'],
            total_spots=sra_summary[uid]['expxml']['Summary']['Statistics']['@total_spots'],
            total_bases=sra_summary[uid]['expxml']['Summary']['Statistics']['@total_bases'],
            total_size=sra_summary[uid]['expxml']['Summary']['Statistics']['@total_size'],
            experiment=sra_summary[uid]['expxml']['Experiment']['@acc'],
            platform=sra_summary[uid]['expxml']['Summary']['Platform']['#text'],
            model=sra_summary[uid]['expxml']['Summary']['Platform']['@instrument_model'],
            tax_id=sra_summary[uid]['expxml']['Organism']['@taxid'],
            sample=sra_summary[uid]['expxml']['Sample']['@acc'],
            library_strategy=sra_summary[uid]['expxml']['Library_descriptor']['LIBRARY_STRATEGY'],
            library_selection=sra_summary[uid]['expxml']['Library_descriptor']['LIBRARY_SELECTION'],
            library_source=sra_summary[uid]['expxml']['Library_descriptor']['LIBRARY_SOURCE'],
            bio_project=sra_summary[uid]['expxml']['Bioproject'],
            bio_sample=sra_summary[uid]['expxml']['Biosample'],
            sra_study=sra_summary[uid]['expxml']['Study']['@acc'],
            sra_id=sra_summary[uid]['expxml']['Submitter']['@acc']
        ))
    return pd.DataFrame(parsed)


def experiment_summary(gse_id):
    """
    Returns an object with the following structure:
    {
        gse_id: str,
        microarray: pd.DataFrame,
        rnaseq: pd.DataFrame if exists
    }
    """
    output = dict(gse_id=gse_id)
    gse_summary = sum_GSE(gse_id)
    output['microarray'] = parse_gse_summary(gse_summary)
    
    # handling rnaseq data
    relevant_srp_ids = find_relevant_srp(gse_summary)
    if relevant_srp_ids:
        sra_ids = entrez_search(db_name=SRA_DB_NAME, term= '+OR+'.join(relevant_srp_ids))
        output['rnaseq'] = parse_sra_summary(entrez_summary(db_name=SRA_DB_NAME, ids=sra_ids)) 
    return output


def save_to_csv(experiment_summary: dict):
    """
    Saves data from an experiment summary to CSV for the purposes of the assignment
    """
    experiment_summary['microarray'].to_csv(MICROARRAY_PATH.format(gse_id=experiment_summary['gse_id']))
    if 'rnaseq' in experiment_summary:
        experiment_summary['rnaseq'].to_csv(RNASEQ_PATH.format(gse_id=experiment_summary['gse_id']))


if __name__ == '__main__':
    save_to_csv(experiment_summary('GSE89408'))
    save_to_csv(experiment_summary('GSE59847'))
    save_to_csv(experiment_summary('GSE40598'))
