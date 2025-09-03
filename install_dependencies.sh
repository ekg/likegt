#!/bin/bash

# LikeGT Dependencies Installation Script
# This script helps install the required dependencies for LikeGT

set -e

echo "======================================"
echo "LikeGT Dependencies Installation"
echo "======================================"
echo

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to check if a command exists
check_command() {
    if command -v $1 &> /dev/null; then
        echo -e "${GREEN}✓${NC} $1 is installed ($(command -v $1))"
        return 0
    else
        echo -e "${RED}✗${NC} $1 is not installed"
        return 1
    fi
}

# Function to install with conda
install_with_conda() {
    echo -e "${YELLOW}Installing $1 with conda...${NC}"
    conda install -c bioconda -y $1 || mamba install -c bioconda -y $1
}

# Function to install with cargo
install_with_cargo() {
    echo -e "${YELLOW}Installing $1 with cargo...${NC}"
    cargo install $1
}

echo "Checking system requirements..."
echo

# Check for package managers
HAS_CONDA=false
HAS_CARGO=false
HAS_APT=false

if command -v conda &> /dev/null || command -v mamba &> /dev/null; then
    HAS_CONDA=true
    echo -e "${GREEN}✓${NC} Conda/Mamba detected"
fi

if command -v cargo &> /dev/null; then
    HAS_CARGO=true
    echo -e "${GREEN}✓${NC} Cargo detected"
else
    echo -e "${YELLOW}!${NC} Cargo not detected. Installing Rust..."
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
    source $HOME/.cargo/env
    HAS_CARGO=true
fi

if command -v apt &> /dev/null; then
    HAS_APT=true
    echo -e "${GREEN}✓${NC} APT detected (Ubuntu/Debian)"
fi

echo
echo "Checking dependencies..."
echo

# Core dependencies
MISSING_DEPS=()

# Check each dependency
echo "Core dependencies:"
check_command "odgi" || MISSING_DEPS+=("odgi")
check_command "minimap2" || MISSING_DEPS+=("minimap2")
check_command "wgsim" || MISSING_DEPS+=("wgsim")
check_command "samtools" || MISSING_DEPS+=("samtools")
check_command "seqtk" || MISSING_DEPS+=("seqtk")
check_command "gfainject" || MISSING_DEPS+=("gfainject")
check_command "gafpack" || MISSING_DEPS+=("gafpack")

echo
echo "Additional dependencies:"
check_command "allwave" || MISSING_DEPS+=("allwave")
check_command "seqwish" || MISSING_DEPS+=("seqwish")
check_command "bwa" || echo "  (optional)"

echo
echo "System tools:"
check_command "zcat" || MISSING_DEPS+=("gzip")
check_command "python3" || MISSING_DEPS+=("python3")

if [ ${#MISSING_DEPS[@]} -eq 0 ]; then
    echo
    echo -e "${GREEN}All dependencies are installed!${NC}"
    echo
    echo "You can now build LikeGT with:"
    echo "  cargo build --release"
    exit 0
fi

echo
echo -e "${YELLOW}Missing dependencies: ${MISSING_DEPS[@]}${NC}"
echo

# Offer installation options
if [ "$HAS_CONDA" = true ]; then
    echo "Would you like to install missing dependencies with conda/mamba? (recommended)"
    echo "This will create a new environment called 'likegt'"
    read -p "Install with conda? [Y/n] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]] || [ -z "$REPLY" ]; then
        # Create conda environment
        echo "Creating conda environment 'likegt'..."
        conda create -n likegt -y python=3.10 || mamba create -n likegt -y python=3.10
        
        echo "Activating environment..."
        eval "$(conda shell.bash hook)"
        conda activate likegt
        
        # Install available tools via conda
        echo "Installing bioinformatics tools..."
        conda install -c bioconda -c conda-forge -y \
            minimap2 \
            samtools \
            seqtk \
            bwa \
            odgi \
            seqwish \
            wgsim || true
        
        # Tools that might need special handling
        echo
        echo "Some tools may need to be installed from source:"
        
        if ! command -v gfainject &> /dev/null; then
            echo "  - gfainject: Install from https://github.com/ekg/gfainject"
        fi
        
        if ! command -v gafpack &> /dev/null; then
            echo "  - gafpack: Install from https://github.com/ekg/gafpack"
        fi
        
        if ! command -v allwave &> /dev/null; then
            echo "  - allwave: Install from https://github.com/ekg/allwave"
        fi
        
        echo
        echo -e "${GREEN}Installation complete!${NC}"
        echo "Remember to activate the environment before using LikeGT:"
        echo "  conda activate likegt"
    fi
elif [ "$HAS_APT" = true ]; then
    echo "Would you like to install available dependencies with apt?"
    read -p "Install with apt? [Y/n] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]] || [ -z "$REPLY" ]; then
        echo "Installing with apt (may require sudo)..."
        sudo apt update
        sudo apt install -y \
            minimap2 \
            samtools \
            seqtk \
            bwa \
            python3 \
            python3-pip \
            build-essential \
            cmake \
            zlib1g-dev
        
        echo
        echo "Some tools need to be installed from source:"
        echo "  - odgi: https://github.com/pangenome/odgi"
        echo "  - seqwish: https://github.com/ekg/seqwish"
        echo "  - gfainject: https://github.com/ekg/gfainject"
        echo "  - gafpack: https://github.com/ekg/gafpack"
        echo "  - allwave: https://github.com/ekg/allwave"
        echo "  - wgsim: https://github.com/lh3/wgsim"
    fi
else
    echo "Manual installation required. Please install the following tools:"
    for dep in "${MISSING_DEPS[@]}"; do
        echo "  - $dep"
    done
    echo
    echo "Installation instructions:"
    echo "  - odgi: https://github.com/pangenome/odgi"
    echo "  - minimap2: https://github.com/lh3/minimap2"
    echo "  - samtools: http://www.htslib.org/"
    echo "  - seqtk: https://github.com/lh3/seqtk"
    echo "  - wgsim: https://github.com/lh3/wgsim"
    echo "  - seqwish: https://github.com/ekg/seqwish"
    echo "  - gfainject: https://github.com/ekg/gfainject"
    echo "  - gafpack: https://github.com/ekg/gafpack"
    echo "  - allwave: https://github.com/ekg/allwave"
fi

echo
echo "After installing dependencies, build LikeGT with:"
echo "  cargo build --release"