import xml.etree.ElementTree as ET
import pandas as pd

# Parse an XML file containing BioSample data, extract specific attributes
# and save the extracted information into a CSV file.
def parse_xml(input, output):
    # Define the columns for the output CSV.
    columns = ['Accession', 'Species', 'Sample_ID', 'Locality', 'Ecotype', 'Tissue']
    extracted_data = []  # Initialize a list to store rows of extracted information.
    
    # Parse the XML file to create an ElementTree structure.
    tree = ET.parse(input)
    root = tree.getroot()  # Access the root element of the XML document.
    
    biosamples = root.findall('BioSample')
    selected_elements = biosamples[:2] + biosamples[41:43] + biosamples[101:103]
    
    # Loop through all 'BioSample' elements in the XML.
    for biosample in selected_elements:
        # Extract elements
        accession = biosample.get('accession', 'Unknown_Accession')
        organism_element = biosample.find('./Description/Organism/OrganismName')
        species = '_'.join(organism_element.text.split())
        
        # Extract the 'ecotype' and 'tissue' attributes with error handling.
        ecotype_element = biosample.find('./Attributes/Attribute[@attribute_name="ecotype"]')
        habitat = ecotype_element.text
                
        # Split 'habitat' into 'locality' and 'ecotype'.
        locality = ' '.join(habitat.split(' ')[:-1])
        ecotype = habitat.split(' ')[-1]

        # Generate and store a unique sample ID
        first_letter = species.split('_')[0][0]
        first_four = species.split('_')[1][:4]
        hab_acronym = 'S' if ecotype == 'Sulfidic' else 'NS'
        sample_id = f"{first_letter}{first_four}{hab_acronym}_{accession}"
        
        tissue_element = biosample.find('./Attributes/Attribute[@attribute_name="tissue"]')
        tissue = tissue_element.text
        
        # Append the extracted data as a dictionary to the list.
        extracted_data.append({
            'Accession': accession,
            'Species': species,
            'Sample_ID': sample_id,
            'Locality': locality,
            'Ecotype': ecotype,
            'Tissue': tissue
        })
    
    # Convert the list of dictionaries into a DataFrame.
    df = pd.DataFrame(extracted_data, columns=columns)
    
    # Write the DataFrame to the specified output CSV file.
    df.to_csv(output, index=False)


# Function to parse sample information from a CSV file. 
# If the file doesn't exist, an empty list is returned.
def get_targets_from_csv(csv_path, suffix):
    try:
        df = pd.read_csv(csv_path)
        # Create paths to target FastQC HTML reports for each sample
        targets = [
            f"data/Genomic_data/{row['Species']}/{row['Accession']}/{row['Accession']}{suffix}"
            for _, row in df.iterrows()
        ]
        return targets
    except FileNotFoundError:
        # Return an empty list if the CSV file is missing
        return []



