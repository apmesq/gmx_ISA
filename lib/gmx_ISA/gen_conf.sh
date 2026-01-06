#!/bin/bash

# Name of the configuration file
config_file="variables.conf"

# Verify if the file exists, if not it creates a new config_file
if [[ ! -f "$config_file" ]]; then
    echo "Creating the configuration file '$config_file'."
    touch "$config_file"
fi

# Function to update or add a new variable on file
update_config() {
    local key=$1
    local value=$2
    if grep -q "^$key=" "$config_file"; then
        sed -i "s/^$key=.*/$key=\"$value\"/" "$config_file"
    else
        echo "$key=\"$value\"" >> "$config_file"
    fi
}

# List of index numbers and variables
declare -A variables=(
    ["Protein"]="What is the number of the term 'Protein' in your index (*.ndx) file?"
    ["Ligand"]="What is the number of the term 'Ligand' in your index (*.ndx) file?"
    ["System"]="What is the number of the term 'System' in your index (*.ndx) file?"
    ["Backbone"]="What is the number of the term 'Backbone' in your index (*.ndx) file?"
    ["Ligand2"]="What is the number of the term 'Ligand2' in your index (*.ndx) file?"
    ["TotalE"]="What is the number for the term 'Total Energy' (it will be followed by '0')?"
    ["CoulombSR"]="What is the number for the term 'Coulomb-SR' (it will be followed by '0')?"
    ["LJSR"]="What is the number for the term 'Lennard-Jones (LJ-SR)' (it will be followed by '0')?"
)

# Iterate over each variable and prompt the user for the value
for key in "${!variables[@]}"; do
    echo "${variables[$key]}"
    read -p "> " value
    # Add " 0" only to TotalE, CoulombSR and LJSR
    if [[ "$key" == "TotalE" || "$key" == "CoulombSR" || "$key" == "LJSR" ]]; then
        value="$value 0"
    fi
    update_config "$key" "$value"
done

echo "The file '$config_file' was written succesfully!"
