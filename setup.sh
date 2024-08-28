#!/bin/bash

# Script to install and run BacPage pipeline

# Function to install mamba
install_mamba() {
    echo "Installing Mamba..."
    curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
    bash Mambaforge-$(uname)-$(uname -m).sh
    echo "Mamba installation complete."
}

# Function to clone the BacPage repository
clone_repo() {
    echo "Cloning BacPage repository..."
    git clone https://github.com/seqafrica/bacpage_nmimr.git
    cd bacpage_nmimr || exit
    git checkout -b split_into_command
    echo "Repository cloned and switched to the development branch."
}

# Function to set up BacPage using conda
setup_conda() {
    echo "Setting up BacPage with Conda..."
    mamba env create -f environment.yaml
    mamba activate bacpage
    pip install .
    echo "BacPage setup complete. Test the installation with 'bacpage -h' and 'bacpage version'."
}

# Function to build Docker image for BacPage
build_docker_image() {
    echo "Building Docker image for BacPage..."
    docker build -t bacpage:latest .
    echo "Docker image built successfully."
}

# Function to run BacPage using Docker
run_docker() {
    echo "Running BacPage with Docker..."
    if [ -z "$1" ]; then
        echo "Please provide the path to your project directory."
        exit 1
    fi
    docker run --rm -it -v "$1":/data bacpage:latest bacpage assemble /data/[your-project-directory-name]
    echo "BacPage run complete. Check the output in your project directory."
}

# Function to update BacPage
update_bacpage() {
    echo "Updating BacPage..."
    git pull
    mamba activate bacpage
    mamba env update -f environment.yaml
    pip install .
    echo "BacPage updated successfully."
}

# Function to update Docker image
update_docker_image() {
    echo "Updating Docker image for BacPage..."
    git pull
    docker build -t bacpage:latest .
    echo "Docker image updated successfully."
}

# Main menu
echo "BacPage Pipeline Setup Script"
echo "Choose an option:"
echo "1. Install Mamba"
echo "2. Clone BacPage repository"
echo "3. Setup BacPage with Conda"
echo "4. Build Docker image for BacPage"
echo "5. Run BacPage with Docker"
echo "6. Update BacPage (Conda)"
echo "7. Update BacPage Docker image"
echo "8. Exit"

read -rp "Enter your choice: " choice

case $choice in
    1) install_mamba ;;
    2) clone_repo ;;
    3) setup_conda ;;
    4) build_docker_image ;;
    5) 
       read -rp "Enter the path to your project directory: " project_path
       run_docker "$project_path" ;;
    6) update_bacpage ;;
    7) update_docker_image ;;
    8) exit 0 ;;
    *) echo "Invalid choice. Exiting." ;;
esac
