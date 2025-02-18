import os
import pandas as pd

# Define your list of compound IDs (cpd numbers)
cpd_list = [
    "cpd09971", "cpd03481", "cpd04533", "cpd21213", 
    "cpd24042", "cpd00363", "cpd10408", "cpd03481",
    "cpd01834", "cpd20521", "cpd05283", 
    "cpd00324", "cpd05284", "cpd01004", "cpd01290",
    "cpd27450", "cpd00150", "cpd00013", "cpd25057",
    "cpd05178", "cpd15644", "cpd05284", "cpd33443",
    "cpd09228", "cpd08349", "cpd17457"
]

# Path to the folder containing SBML files
sbml_folder = "./16-carveme"

# Dictionary to store results
cpd_results = {cpd: [] for cpd in cpd_list}

# Loop through all SBML (XML) files in the folder
for filename in os.listdir(sbml_folder):
    if filename.endswith(".xml"):  # Ensure it's an SBML file
        sbml_path = os.path.join(sbml_folder, filename)

        try:
            # Read file as plain text
            with open(sbml_path, "r", encoding="utf-8") as f:
                content = f.read()

                # Search for each cpd number
                for cpd in cpd_list:
                    if cpd in content:  # If found, add to results
                        cpd_results[cpd].append(filename)

        except Exception as e:
            print(f"Error processing {filename}: {e}")

# Convert results to a Pandas DataFrame
df = pd.DataFrame([(cpd, ", ".join(set(files))) for cpd, files in cpd_results.items()], columns=["Compound", "SBML Files Found In"])

# Define file paths for metadata
metadata_file_1 = "./01-metadata/id_to_minionbarcode.txt"  # Contains: Genus, internalID, barcode
metadata_file_2 = "./01-metadata/isolateinfo.txt"  # Contains: internalID, species

# Load metadata files
df_meta1 = pd.read_csv(metadata_file_1, sep="\t")  # Adjust separator if needed
df_meta1["barcode"] = df_meta1["barcode"].astype(str).str.zfill(2)
df_meta2 = pd.read_csv(metadata_file_2, sep="\t")  

# Merge metadata to map barcode -> species name
df_metadata = df_meta1.merge(df_meta2, on="shortIDreplacement", how="left")  # Merge on internalID


# Create a dictionary to map barcodes to species names
barcode_to_species = dict(zip(df_metadata["barcode"], df_metadata["SpeciesWGS"]))

# Fix SBML filenames → Extract barcode number before mapping
def map_barcodes_to_species(file_list):
    species_list = set()
    
    for filename in file_list.split(", "):
        barcode_number = filename.replace("barcode", "").replace(".xml", "").zfill(2)  # Extract numeric part
        species = barcode_to_species.get(barcode_number, filename)  # Map to species name
        species_list.add(species)

    return ", ".join(species_list)

# Apply transformation
df["SBML Files Found In"] = df["SBML Files Found In"].apply(map_barcodes_to_species)


# Save updated results
updated_output_file = "./16-carveme/cpd_search_VOCs.csv"
df.to_csv(updated_output_file, index=False)

