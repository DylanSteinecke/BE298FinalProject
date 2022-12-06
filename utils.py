import pandas as pd
import matplotlib.pyplot as plt, numpy as np, random
import xml.etree.ElementTree as ET
import json, requests as req, os
from IPython.core.display import display, HTML
display(HTML("<style>.container { width:90% !important; }</style>"))

def switch_dictset_to_dictlist(the_dict):
    '''
    FUNCTION:
    - Make a new dictionary with values as lists 
      instead of values as sets
      
    PARAMS:
    - the_dict: The initial dict with values of sets
    '''
    
    dictlist = dict()
    
    for k in the_dict.copy():
        dictlist[k] = list(the_dict[k])
        
    return dictlist


def map_drug_to_go_term_proteins(relation, drug_to_protein_dicts, 
                                 go_proteins, go_term):
    '''
    FUNCTION:
    - Map drugs to their go proteins
    
    PARAMS:
    - relation: 'target', 'carrier', 'enzyme', or 'transporter'
    - go_proteins: proteins associated with the GO term
    - go_term: the name of the GO term, for writing to file
    '''
    
    drug_to_go_proteins = dict()

    # Drug
    for drug, proteins in drug_to_protein_dicts[relation].items():
        drug_to_go_proteins[drug] = set()

        # Protein
        for protein in proteins:
            # GO term-related protein
            if protein in go_proteins:
                drug_to_go_proteins[drug].add(protein)

    drug_to_go_proteins = switch_dictset_to_dictlist(drug_to_go_proteins)
    outpath = f'data/drug_to_{go_term}_{relation}s.json'
    json.dump(drug_to_go_proteins, open(outpath,'w'))
    
    return drug_to_go_proteins



def get_drugs_proportion_go_proteins(drug_to_protein, drug_to_go_protein):
    '''
    FUNCTION:
    - Takes a drug-[relation]-protein dict and finds the proportion of
      proteins associated with your GO term
      
    PARAMS:
    '''
    drug_to_proportion_go_proteins = dict() 

    for drug, go_proteins in drug_to_go_protein.items():
        num_go_proteins = len(go_proteins)
        num_all_proteins = len(drug_to_protein[drug])
        drug_to_proportion_go_proteins[drug] = num_go_proteins/num_all_proteins
        
    drug_to_proportion_go_proteins = sorted(drug_to_proportion_go_proteins.items(), key=lambda x:x[1], reverse=True)
        
    return drug_to_proportion_go_proteins


def get_protein_name_to_ids():
    protein_name_to_id_req = req.get('https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cprotein_name&format=tsv&query=%28%2A%29%20AND%20%28reviewed%3Atrue%29%20AND%20%28model_organism%3A9606%29')
    n_split = protein_name_to_id_req.text.split('\n')
    list_of_lists = [row.split('\t') for row in n_split]
    protein_name_df = pd.DataFrame(list_of_lists, columns=list_of_lists[0])
    protein_name_df.drop(0, inplace=True)
    protein_name_df.drop(len(protein_name_df), inplace=True)

    all_protein_names = list()
    protein_id_to_name, protein_name_to_id = dict(), dict()

    for protein_id, protein_names in zip(protein_name_df['Entry'], protein_name_df['Protein names']):
        try: protein_name = protein_names.split('(')[0]
        except: print(protein_names)
        all_protein_names.append(protein_name)
        protein_id_to_name[protein_id] = protein_name
        protein_name_to_id[protein_name] = protein_id

    return protein_name_to_id, protein_id_to_name


def get_gene_disease_association_df(download=False):
    if download:
        os.terminal('wget -N -P data/ https://www.disgenet.org/static/disgenet_ap1/files/downloads/all_gene_disease_associations.tsv.gz')
        os.terminal('gunzip data/all_gene_disease_associations.tsv.gz')
    gda_df = pd.read_table('data/all_gene_disease_associations.tsv')
    return gda_df
    
def get_disease_id_to_name(gda_df):
    disease_id_to_name, disease_name_to_id = dict(), dict()
    for disease_id, disease_name in zip(gda_df['diseaseId'], gda_df['diseaseName']):
        disease_id_to_name[disease_id] = disease_name
        disease_name_to_id[disease_name] = disease_id
    return disease_id_to_name, disease_name_to_id


def get_protein_disease_associations(gda_df):
    # Gene/Protein -[associated with]- Disease
    kgdr_path = '../Drug Repurposing Knowledge Graph/'
    protein2gene = json.load(open(kgdr_path+'output/protein2gene/all_uniprot2entrez.json'))
    gene2protein = json.load(open(kgdr_path+'output/protein2gene/all_entrez2uniprot.json'))
    gene2disease, disease2gene, protein2disease, disease2protein = dict(), dict(), dict(), dict()
    for gene, disease in zip(gda_df['geneId'], gda_df['diseaseId']):
        gene2disease.setdefault(gene, set()).add(disease)
        disease2gene.setdefault(disease, set()).add(gene)
        try:
            proteins = gene2protein[str(gene)]
            for protein in proteins:
                protein2disease.setdefault(protein, set()).add(disease)
                disease2protein.setdefault(disease, set()).add(protein)
        except:
            continue
    #gene2disease = dict(sorted(gene2disease.items(), key=lambda x:x[1], reverse=True))
    #disease2gene = dict(sorted(disease2gene.items(), key=lambda x:x[1], reverse=True))

    return protein2disease, disease2protein, gene2disease, disease2gene



def get_disease_to_protein_dicts(disease2protein, human_go_proteins, disease_id_to_name, protein_id_to_name):
    disease_to_go_protein, go_protein_to_disease = dict(), dict()
    disease_name_to_go_protein_name, go_protein_name_to_disease_name = dict(), dict()

    for disease_id, proteins in disease2protein.items():
        for protein_id in proteins:
            if protein_id in human_go_proteins:

                # Disease Name, Protein Name
                disease_name = disease_id_to_name[disease_id]
                try: protein_name = protein_id_to_name[protein_id]
                except: protein_name = protein_id

                # Disease -[associated with]- Protein
                disease_to_go_protein.setdefault(disease_id, set()).add(protein_id)
                go_protein_to_disease.setdefault(protein_id, set()).add(disease_id)
                disease_name_to_go_protein_name.setdefault(protein_name, set()).add(disease_name)
                go_protein_name_to_disease_name.setdefault(disease_name, set()).add(protein_name)
 
    disease_to_go_protein = switch_dictset_to_dictlist(disease_to_go_protein)
    disease_name_to_go_protein_name = switch_dictset_to_dictlist(disease_name_to_go_protein_name)
    go_protein_to_disease = switch_dictset_to_dictlist(go_protein_to_disease)
    disease_name_to_go_protein_name = switch_dictset_to_dictlist(disease_name_to_go_protein_name)

    return disease_to_go_protein, go_protein_to_disease, disease_name_to_go_protein_name, go_protein_name_to_disease_name



def use_high_quality_sources_gda(gda_df):
    ''' gda_df = gene disease association dataframe'''
    high_quality_sources = ['GWASDB','GWASCAT','CTD_human','ORPHANET','CLINVAR','HPO','CLINGEN','GENOMICS ENGLAND','PSYGENET','UNIPROT']
    high_quality_rows = list()
    for row_num, source_entry in enumerate(gda_df['source']):
        if source_entry in high_quality_sources:
            high_quality_rows.append(True)
        else:
            high_quality_rows.append(False)
    df = gda_df[high_quality_rows]
    return df



def get_mesh_name_id_tree_mappings(root, download=False):
    
    if download:
        os.terminal('wget -N -P data/ https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/desc2022.xml')

    name_to_id, id_to_name, id_to_tree, tree_to_id = dict(), dict(), dict(), dict()
    all_tree_numbers = list()

    for ele in root:
        try:
            # MeSH Tree Number
            tree_numbers = ele.find('TreeNumberList').findall('TreeNumber')
            for tree_number in tree_numbers:
                # Disease tree number
                if True:#tree_number.text.startswith(('C','F03')):

                    # Tree
                    tree_number = tree_number.text
                    all_tree_numbers.append(tree_number)

                    # ID to Tree
                    try:
                        ID = ele.find('DescriptorUI').text
                        id_to_tree.setdefault(ID,set()).add(tree_number)
                        tree_to_id.setdefault(tree_number,set()).add(ID)
                    except:
                        pass

                    # ID to Name
                    try:
                        name = ele.find('DescriptorName').find('String').text
                        name_to_id.setdefault(name,set()).add(ID)
                        id_to_name.setdefault(ID,set()).add(name)
                    except:
                        pass
        except:
            continue        


    all_tree_numbers = sorted(all_tree_numbers)
    tree_to_id = dict(sorted(tree_to_id.items()))

    name_to_id = switch_dictset_to_dictlist(name_to_id)
    id_to_name = switch_dictset_to_dictlist(id_to_name)
    tree_to_id = switch_dictset_to_dictlist(tree_to_id)
    id_to_tree = switch_dictset_to_dictlist(id_to_tree)

    return name_to_id, id_to_name, tree_to_id, id_to_tree, all_tree_numbers


def get_mesh_tree_hierarchy(all_tree_numbers):
    tree_to_tree = dict()

    # Tree Number
    for tree_num in all_tree_numbers:
        if '.' in tree_num:

            # Parent of Tree Number
            parent = ''
            for num in tree_num.split('.')[:len(tree_num.split('.'))-1]:
                parent += num+'.'
            parent = parent.strip('.')

            # Tree Number -[subclass of]-> Tree Number
            tree_to_tree[tree_num] = parent
            
    return tree_to_tree


def get_umls_disease_category(disease_category, umls_to_mesh_id):
    '''
    FUNCTION:
    - Figure out which UMLS disease IDs belong to a disease category
      of your choice. This maps UMLS to MeSH tree numbers to find 
      the disease category. 
    PARAMS:
    - disease_category: The prefix of the MeSH tree number category
    '''
    umls_disease_in_category = set()
    for umls_disease in diseases_with_most_mitochondrial_proteins:   
        try: mesh_disease_ids = umls_to_mesh_id[umls_disease]
        except: continue
        for mesh_disease_id in mesh_disease_ids:
            try: mesh_disease_trees = mesh_id_to_tree[mesh_disease_id]
            except: continue
            for mesh_disease_tree in mesh_disease_trees:
                if mesh_disease_tree.startswith(disease_category):
                    umls_disease_in_category.add(umls_disease)
    return umls_disease_in_category

def get_umls_mesh_category(mesh_category, umls_ids, umls_to_mesh_id, mesh_id_to_tree):
    '''
    FUNCTION:
    - Figure out which UMLS disease IDs belong to a disease category
      of your choice. This maps UMLS to MeSH tree numbers to find 
      the disease category. 
    PARAMS:
    - disease_category: The prefix of the MeSH tree number category
    '''
    umls_in_category = set()
    # UMLS
    for umls_id in umls_ids:   
        # UMLS-MeSH
        try: mesh_ids = umls_to_mesh_id[umls_id]
        except: continue
        for mesh_id in mesh_ids:
            # MeSH - MeSH Tree
            try: mesh_trees = mesh_id_to_tree[mesh_id]
            except: continue
            for mesh_tree in mesh_trees:
                # Correct MeSH Tree
                if mesh_tree.startswith(mesh_category):
                    umls_in_category.add(umls_id)
    return umls_in_category

                    
    
def get_drug_treats_disease(db_to_mesh_df, mesh_to_umls, umls_your_diseases):
    drug_treats_disease, disease_treated_drug = dict(), dict()
    for drug, mesh_disease in zip(db_to_mesh_df['Compound (DrugBank)'], db_to_mesh_df['Disease (MeSH)']):
        drug = drug.split(':')[1]
        mesh_disease = mesh_disease.split(':')[1]

        try: umls_diseases = mesh_to_umls[mesh_disease]
        except: continue
        for umls_disease in umls_diseases:
            drug_treats_disease.setdefault(drug, set()).add(umls_disease)
            disease_treated_drug.setdefault(umls_disease, set()).add(drug)
    drug_treats_disease = switch_dictset_to_dictlist(drug_treats_disease)
    disease_treated_drug = switch_dictset_to_dictlist(disease_treated_drug)

    your_treated_diseases = list(set(disease_treated_drug).intersection(umls_your_diseases))
    
    return drug_treats_disease, disease_treated_drug, your_treated_diseases




def get_drugs_predicted_phenotypes(result_table_dir, bh_cutoff):
    '''
    FUNCTION:
    - Get each drug's predicted phenotypes
    PARAMS:
    - bh_cutoff: Benjamini-Hochberg cutoff for significance
    '''
    predicted_drug_to_disease = dict()
    for file_name in os.listdir(result_table_dir):
        drug_id = file_name
        result_link = f'PathFX/results/all_network_results/{drug_id}/{drug_id}_merged_neighborhood__assoc_table_.txt'
        try: drug_disease_prediction_df = pd.read_table(result_link)
        except: continue
        if len(drug_disease_prediction_df) > 1:
            df = drug_disease_prediction_df.copy()
            df = df[df['Benjamini-Hochberg']<bh_cutoff]
            disease_ids = list(df['cui'])
            for disease_id in disease_ids:
                predicted_drug_to_disease.setdefault(drug_id, set()).add(disease_id)
    
    predicted_drug_to_disease = switch_dictset_to_dictlist(predicted_drug_to_disease)
    return predicted_drug_to_disease



def get_pathways_genes_and_proteins(correct_drug_to_disease, bh_cutoff):
    '''
    FUNCTION:
    - Get the genes/proteins involved in the PPI network path between
      the drug and the disease for the correctly identified drug-disease pairs.
    PARAMS:
    - correct_drug_to_disease: The drugs that are correctly identified in drug-disease
      pairs generated by PathFX
    - bh_cutoff: Benjamini-Hochberg cutoff for significance
    '''

    drug_to_disease_to_genes = dict()
    drug_to_disease_to_proteins = dict()

    for drug_id in correct_drug_to_disease:
        result_link = f'PathFX/results/all_network_results/{drug_id}/{drug_id}_merged_neighborhood__assoc_table_.txt'
        drug_disease_prediction_df = pd.read_table(result_link)
        if len(drug_disease_prediction_df) > 1:
            df = drug_disease_prediction_df.copy()
            df = df[df['Benjamini-Hochberg']<bh_cutoff]

            # Drug-Disease-Gene Name
            disease_to_genes = dict()
            disease_to_proteins = dict()
            for disease_id, genes in zip(df['cui'], df['genes']):
                genes = genes.split(',')
                disease_to_genes[disease_id] = genes

                # Drug-Disease-Protein ID
                proteins = list()
                for gene_name in genes:
                    try: 
                        protein_ids = [protein_id for protein_id in gene_name_to_protein_id[gene_name]]   
                        proteins += protein_ids
                    except: 
                        print(f'Couldn\'t map {drug}\'s gene {gene_name} for {disease_id}')
                disease_to_proteins[disease_id] = proteins

            drug_to_disease_to_genes[drug_id] = disease_to_genes
            drug_to_disease_to_proteins[drug_id] = disease_to_proteins
            
    return drug_to_disease_to_genes, drug_to_disease_to_proteins



def map_protein_ids_to_gene_names():
    # Download table of mappings
    r = req.get('https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cgene_names%2Cgene_primary%2Cgene_synonym&format=tsv&query=%28%2A%29%20AND%20%28model_organism%3A9606%29%20AND%20%28reviewed%3Atrue%29')
    n_split = r.text.split('\n')
    list_of_lists = [row.split('\t') for row in n_split]
    gene_name_df = pd.DataFrame(list_of_lists, columns=list_of_lists[0])
    gene_name_df.drop(0, inplace=True)
    gene_name_df.drop(len(gene_name_df), inplace=True)

    # Extract mappings, load into dictionaries
    protein_id_to_gene_name, gene_name_to_protein_id = dict(), dict()
    for protein_id, gene_names in zip(gene_name_df['Entry'], gene_name_df['Gene Names']):
        gene_names = gene_names.split(' ')
        for gene_name in gene_names:
            protein_id_to_gene_name.setdefault(protein_id, set()).add(gene_name)
            gene_name_to_protein_id.setdefault(gene_name, set()).add(protein_id)
    gene_name_to_protein_id = switch_dictset_to_dictlist(gene_name_to_protein_id)
    protein_id_to_gene_name = switch_dictset_to_dictlist(protein_id_to_gene_name)
    
    return gene_name_to_protein_id, protein_id_to_gene_name



def get_human_proteins():
    try:
        human_protein_ids = set(json.load(open('data/human_protein_ids.json')))
    except:
        human_protein_ids = set(req.get('https://rest.uniprot.org/uniprotkb/stream?format=list&query=%28%28proteome%3AUP000005640%29%29%20AND%20%28reviewed%3Atrue%29').text.split('\n'))
        json.dump(list(human_protein_ids), open('data/human_protein_ids.json','w'))
    return human_protein_ids    


def get_all_drugs_proteins(drug_to_protein_dict):
    all_drugs_proteins = list()
    for proteins in drug_to_protein_dict.values():
        all_drugs_proteins += proteins
        
    unique_all_drugs_proteins = set(all_drugs_proteins)
    return unique_all_drugs_proteins, all_drugs_proteins


def get_mitochondrial_drugs_proteins(all_drugs_proteins, mitochondrial_proteins):
    mitochondrial_drugs_proteins, unique_mitochondrial_drugs_proteins = list(), set()

    for drugs_proteins in all_drugs_proteins:
        if drugs_proteins in mitochondrial_proteins:
            unique_mitochondrial_drugs_proteins.add(drugs_proteins)
            mitochondrial_drugs_proteins.append(drugs_proteins)
    return unique_mitochondrial_drugs_proteins, mitochondrial_drugs_proteins  


def sort_dict_by_values_length(the_dict):
    return sorted(the_dict, key = lambda k:len(the_dict[k]), reverse=True)

def sort_dict_by_values(the_dict):
    return sorted(the_dict, key = lambda k:the_dict[k], reverse=True)

def get_drug_treats_your_disease(drug_treats_disease, treated_your_diseases):
    drug_treats_your_diseases = dict()
    your_disease_treated_by_drugs = dict()
    for drug, diseases in drug_treats_disease.items():
        treats_your_diseases = list(set(diseases).intersection(treated_your_diseases))
        if len(treats_your_diseases) > 0:
            drug_treats_your_diseases[drug] = treats_your_diseases
            for your_disease in treats_your_diseases:
                your_disease_treated_by_drugs.setdefault(your_disease, list()).append(drug)
    return your_disease_treated_by_drugs, drug_treats_your_diseases


def get_correctly_predicted_drug_to_disease(predicted_drug_to_disease, drug_treats_disease):
    correct_drug_to_disease = dict()
    for drug, diseases_actual in drug_treats_disease.items():
        try: predicted_diseases = predicted_drug_to_disease[drug]
        except: continue
        correct_predict_diseases = set(diseases_actual).intersection(set(predicted_diseases))
        correct_drug_to_disease[drug] = correct_predict_diseases
    correct_drug_to_disease = switch_dictset_to_dictlist(correct_drug_to_disease)
    
    return correct_drug_to_disease


def get_your_diseases_drugs_proteins(drug_treats_your_disease, drug_to_protein_dict):
    your_disease_drugs_proteins = list()
    for drug in drug_treats_your_disease:
        try:  
            drugs_proteins = drug_to_protein_dict[drug]
            your_disease_drugs_proteins += drugs_proteins
        except: 
            continue
            
    unique_your_disease_drugs_proteins = set(your_disease_drugs_proteins)
    return unique_your_disease_drugs_proteins, your_disease_drugs_proteins


def get_all_pathway_proteins(drug_to_disease_to_proteins):
    '''
    FUNCTION:
    - Get all the proteins involved in the drug pathways
    '''
    all_pathway_proteins = list()
    for drug, disease_to_proteins in drug_to_disease_to_proteins.items():
        for disease, proteins in disease_to_proteins.items():
            all_pathway_proteins += proteins
    all_pathway_unique_proteins = set(all_pathway_proteins)
    
    return all_pathway_unique_proteins, all_pathway_proteins


def get_pathway_mitochondrial_proteins(all_pathway_proteins, mitochondrial_proteins):
    pathway_mitochondrial_proteins = list()
    for pathway_protein in all_pathway_proteins:
        if pathway_protein in mitochondrial_proteins:
            pathway_mitochondrial_proteins.append(pathway_protein)

    pathway_unique_mitochondrial_proteins = set(pathway_mitochondrial_proteins)
    return pathway_unique_mitochondrial_proteins, pathway_mitochondrial_proteins


def display_go_related_proportions(unique_proteins, proteins, 
                                   unique_go_related_proteins, go_related_proteins,
                                   rel, go_term, drug_type, display_it=True):
    #unique_proteins, proteins = get_drugs_proteins(drug_to_protein_dicts[rel[:-1]])
    #unique_go_related_proteins, go_related_proteins = get_go_related_drugs_proteins(proteins, go_related_proteins)
        
    go_related_proportion = round(len(go_related_proteins)/len(proteins),3)*100
    unique_go_related_proportion = round(len(unique_go_related_proteins)/len(unique_proteins),3)*100
        
    if display_it:
        print(f'{len(unique_proteins)} unique {rel} of {len(proteins)} {rel} ')
        print(f'{len(go_related_proteins)} {go_term} {drug_type} {rel}'+\
              f' of all {len(proteins)} {drug_type} {rel}'+\
              f' ({go_related_proportion}%))')
        print(f'{len(unique_go_related_proteins)} unique {go_term} {drug_type} {rel}'+\
          f' of all {len(unique_proteins)} unique {drug_type} {rel}'+\
          f' ({unique_go_related_proportion}%))\n')
    
    
    return go_related_proportion, unique_go_related_proportion


def display_go_random_proportions_sim(unique_proteins, proteins, 
                                   unique_go_random_proteins, go_random_proteins,
                                   go_random_proportion, unique_go_random_proportion,
                                   rel, go_term, drug_type):
    '''
    - unique_proteins: set of relevant proteins (e.g., targets)
    - proteins: list of relevant proteins (e.g., targets)
    - unique_go_random_proportion_proteins: set of "GO-related" proteins 
    - go_random_proportion_proteins: list of random "GO-related" proteins 
    - go_term: GO term name to be displayed (e.g., "Mitochondria")
    '''

    print(f'{len(unique_proteins)} unique {rel}s of {len(proteins)} {rel} ')
    print(f'{len(go_random_proteins)} {go_term} {drug_type} {rel}'+\
          f' of all {len(proteins)} {drug_type} {rel}'+\
          f' ({go_random_proportion}%))')
    print(f'{len(unique_go_random_proteins)} unique {go_term} {drug_type} {rel}s'+\
      f' of all {len(unique_proteins)} unique {drug_type} {rel}s'+\
      f' ({unique_go_random_proportion}%))\n')
    
    
def simulate_go_term_proteins_in_drug_proteins(proteome, original_go_term_proteins, simulations, rel,
                                               go_term, drug_type, specific_drugs = [], all_drugs = True):
    # Note: to test for significance, we should return the distribution, not just the proportions
    '''
    - proteome: list of all proteins in a species (e.g., human proteins)
    - original_go_term_proteins: list of all proteins related to you GO term of interest (e.g., mitochondrial proteins)
    - simulations: number of times to simulate the GO-term-proteins in drug-proteins
    - rel: drug to protein relationship (e.g., drugproteins)
    - go_term: name of the GO term to display
    - drug_type: name of the drug type to display (e.g., cardiovacsular disease drug)
    - specific_drugs: the drugs you care about, you want their proteins (e.g., their targets)
    - all_drugs: True or False if you want all drugs or just some drugs
    '''
    # Simulate a number of times
    sim_targ_cvd_prop, sim_targ_cvd_uniq_prop = 0, 0
    for sim in range(simulations):
        # Random proteins
        fake_mito_proteins = random.sample(list(proteome), len(original_go_term_proteins))
        
        # All drugs or some drugs?
        if all_drugs:
            unique_drugproteins, drugproteins = get_all_drugs_proteins(drug_to_protein_dicts[rel])
        else:
            unique_drugproteins, drugproteins = get_your_diseases_drugs_proteins(specific_drugs, drug_to_protein_dicts[rel])
        
        # Get proportions of GO-term-proteins in drug-proteins
        unique_mitochondrial_drugproteins, mitochondrial_drugproteins = get_mitochondrial_drugs_proteins(drugproteins, fake_mito_proteins)
        temp1, temp2 = display_go_related_proportions(unique_drugproteins, drugproteins, 
                                           unique_mitochondrial_drugproteins, mitochondrial_drugproteins,
                                           rel = rel, go_term = go_term, drug_type = drug_type, display_it = False)
        sim_targ_cvd_prop += temp1
        sim_targ_cvd_uniq_prop += temp2
    
    # Display final average of simulations
    sim_targ_cvd_prop /= simulations
    sim_targ_cvd_uniq_prop /= simulations
    display_go_random_proportions_sim(unique_drugproteins, drugproteins, 
                                   unique_mitochondrial_drugproteins, mitochondrial_drugproteins,
                                   sim_targ_cvd_prop, sim_targ_cvd_uniq_prop,
                                   rel = rel, go_term = go_term, drug_type = drug_type)
    
    
def display_display_disease_to_go_term_proteins(diseases_to_go_related_protein, diseases_to_protein, 
                                               disease_name, go_term):
    diseases_to_number_go_related_proteins = dict()
    cvd_to_proportion_go_related_proteins = dict()
    for cvd, go_related_proteins in diseases_to_go_related_protein.items():
        proteins = diseases_to_protein[cvd]
        diseases_to_number_go_related_proteins[cvd] = len(go_related_proteins)
        cvd_to_proportion_go_related_proteins[cvd] = len(go_related_proteins)/len(proteins)

    # Number of GO Term-Related proteins per disease
    plt.title(f'Number of {go_term} Proteins in Number of {disease_name}s')
    plt.xlabel(f'Number of  {go_term} Proteins')
    plt.ylabel(f'Number of {disease_name}')
    plt.hist(diseases_to_number_go_related_proteins.values(), bins=40);
    
    return cvd_to_proportion_go_related_proteins, diseases_to_number_go_related_proteins

def display_display_disease_to_go_term_proteins_proportion(disease_to_proportion_go_related_proteins, 
                                                           disease_name, go_term):
    # Proportion of GO Term-Related proteins per disease
    plt.title(f'Proportion of {go_term} Proteins in Number of {disease_name}s')
    plt.xlabel(f'Proportion of {go_term} Proteins')
    plt.ylabel(f'Proportion of {disease_name}s')
    plt.hist(disease_to_proportion_go_related_proteins.values(), bins=20);
    
    
    
def simulate_disease_to_go_term_related_protein(diseases_to_go_protein, go_related_proteins,
                                                all_protein_ids, simulations):
    '''
    FUNCTION:
    - Simulate what the disease_to_protein dictionaries could be if the 
      proteins were chosen at random. Then look and see if the subset of
      proteins (e.g., mitochondrial proteins) is different from the 
      actual (i.e., non-simulated) disease_to_protein dictionary
    '''
    simulated_disease_has_mito = 0
    
    for sim in range(simulations):
        diseases_to_random_proteins = dict()
        for disease, proteins in diseases_to_go_protein.items():
            num_proteins = len(proteins)
            random_proteins = set(random.sample(list(all_protein_ids), num_proteins))
            mito_in_random = random_proteins.intersection(go_related_proteins)
            if len(mito_in_random) > 0:
                simulated_disease_has_mito += 1
                
    simulated_disease_has_mito = simulated_disease_has_mito/simulations
                
    return simulated_disease_has_mito


def check_if_disease_to_proteins_is_enriched_for_go_term(cvd_to_protein, go_term):
    '''
    FUNCTION:
    - Check if proteins are enriched for GO term in each disease
    '''
    disease_enriched_for_go_term = list()
    for disease_count, (disease, proteins) in enumerate(disease_to_protein.items()):
        #print(disease_count, end='\r')
        gp = GProfiler(return_dataframe=True)
        df = gp.profile(organism='hsapiens', query=list(proteins))
        if go_term in list(df['name']):
            disease_enriched_for_go_term.append(disease)
    return disease_enriched_for_go_term


def get_drug_disease_enriched_for_go_term(drug_to_disease_to_proteins, go_term = 'mitochondria'):
    drug_disease_enriched_for_go_term = dict()
    for drug, disease_to_proteins in drug_to_disease_to_proteins.items():
        for disease, proteins in disease_to_proteins.items():
            gp = GProfiler(return_dataframe=True)
            df = gp.profile(organism='hsapiens', query=list(proteins))
            if go_term in list(df['name']):
                drug_disease_enriched_for_go_term.setdefault(drug, list()).append(disease)
                
    return drug_disease_enriched_for_go_term


def map_disease_to_drug_to_proteins(drug_to_disease_to_proteins):
    disease_to_drug_to_proteins = dict()
    for drug, disease_to_proteins in drug_to_disease_to_proteins.items():
        for disease, proteins in disease_to_proteins.items() :
            disease_to_drug_to_proteins.setdefault(disease, dict())
            disease_to_drug_to_proteins[disease][drug] = proteins
    return disease_to_drug_to_proteins


def map_proteins_expressed_in_tissue(gene_to_anatomy, gene2protein):
    protein_to_anatomy = dict()
    for gene, anatomies in gene_to_anatomy.items():
        try: proteins = gene2protein[gene]
        except: continue
        for protein in proteins:
            for anatomy in anatomies:
                protein_to_anatomy.setdefault(protein, set()).add(anatomy)
    protein_to_anatomy = switch_dictset_to_dictlist(protein_to_anatomy)

    return protein_to_anatomy
