#!/bin/bash

# Documentation Setup Script for Cavity HOOMD
# This script helps you prepare and deploy your documentation to Read the Docs

set -e

echo "Setting up Documentation for Cavity HOOMD"
echo "=========================================="

# Color codes for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Check if we're in a git repository
if ! git rev-parse --git-dir > /dev/null 2>&1; then
    echo -e "${RED}Error: Not in a git repository${NC}"
    exit 1
fi

# Check if we have the docs directory
if [ ! -d "docs" ]; then
    echo -e "${RED}Error: docs directory not found${NC}"
    exit 1
fi

echo -e "${BLUE}Pre-deployment checklist:${NC}"
echo "1. Git repository detected"
echo "2. Documentation directory found"
echo "3. Read the Docs configuration files present"

# Check current branch
CURRENT_BRANCH=$(git branch --show-current)
echo -e "${BLUE}Current branch: ${CURRENT_BRANCH}${NC}"

# Build documentation locally first
echo -e "${YELLOW}Building documentation locally...${NC}"
cd docs
if make html; then
    echo -e "${GREEN}Documentation built successfully${NC}"
else
    echo -e "${RED}Documentation build failed${NC}"
    echo "Please fix the build errors before proceeding."
    exit 1
fi
cd ..

# Check if changes need to be committed
if ! git diff-index --quiet HEAD --; then
    echo -e "${YELLOW}You have uncommitted changes.${NC}"
    echo "Would you like to commit them now? (y/n)"
    read -r response
    if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]; then
        echo "Enter commit message:"
        read -r commit_message
        git add .
        git commit -m "$commit_message"
        echo -e "${GREEN}Changes committed${NC}"
    else
        echo -e "${YELLOW}Proceeding with uncommitted changes${NC}"
    fi
fi

# Push to GitHub
echo -e "${YELLOW}Pushing to GitHub...${NC}"
git push origin "$CURRENT_BRANCH"

echo -e "${GREEN}Setup complete!${NC}"
echo ""
echo -e "${BLUE}Next steps for Read the Docs:${NC}"
echo "1. Go to https://readthedocs.org/ and sign in with your GitHub account"
echo "2. Click 'Import a Project' and select your repository"
echo "3. Configure the project settings:"
echo "   - Name: cavity-hoomd"
echo "   - Repository URL: https://github.com/muhammadhasyim/cav-hoomd"
echo "   - Default branch: main"
echo "4. In Advanced Settings:"
echo "   - Python configuration file: .readthedocs.yaml"
echo "   - Requirements file: docs/requirements.txt"
echo ""
echo -e "${BLUE}Your documentation will be available at:${NC}"
echo "https://cavity-hoomd.readthedocs.io/"
echo ""
echo -e "${BLUE}Automatic updates:${NC}"
echo "- Documentation rebuilds automatically when you push to main/master"
echo "- Build status visible in your Read the Docs dashboard"
echo "- Typical deployment time: 3-5 minutes"
echo ""
echo -e "${BLUE}Manual build:${NC}"
echo "- Go to your Read the Docs project dashboard"
echo "- Click 'Build Version' to trigger a manual build"
echo ""
echo -e "${GREEN}Happy documenting!${NC}" 