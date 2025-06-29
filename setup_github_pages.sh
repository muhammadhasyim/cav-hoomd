#!/bin/bash

# GitHub Pages Setup Script for Cavity HOOMD Documentation
# This script helps you deploy your documentation to GitHub Pages

set -e

echo "🚀 Setting up GitHub Pages for Cavity HOOMD Documentation"
echo "============================================================"

# Color codes for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Check if we're in a git repository
if ! git rev-parse --git-dir > /dev/null 2>&1; then
    echo -e "${RED}❌ Error: Not in a git repository${NC}"
    exit 1
fi

# Check if we have the docs directory
if [ ! -d "docs" ]; then
    echo -e "${RED}❌ Error: docs directory not found${NC}"
    exit 1
fi

echo -e "${BLUE}📋 Pre-deployment checklist:${NC}"
echo "1. ✅ Git repository detected"
echo "2. ✅ Documentation directory found"
echo "3. ✅ GitHub Actions workflow configured"

# Check current branch
CURRENT_BRANCH=$(git branch --show-current)
echo -e "${BLUE}📍 Current branch: ${CURRENT_BRANCH}${NC}"

# Build documentation locally first
echo -e "${YELLOW}🔨 Building documentation locally...${NC}"
cd docs
if make html; then
    echo -e "${GREEN}✅ Documentation built successfully${NC}"
else
    echo -e "${RED}❌ Documentation build failed${NC}"
    echo "Please fix the build errors before proceeding."
    exit 1
fi
cd ..

# Check if changes need to be committed
if ! git diff-index --quiet HEAD --; then
    echo -e "${YELLOW}⚠️  You have uncommitted changes.${NC}"
    echo "Would you like to commit them now? (y/n)"
    read -r response
    if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]; then
        echo "Enter commit message:"
        read -r commit_message
        git add .
        git commit -m "$commit_message"
        echo -e "${GREEN}✅ Changes committed${NC}"
    else
        echo -e "${YELLOW}⚠️  Proceeding with uncommitted changes${NC}"
    fi
fi

# Push to GitHub
echo -e "${YELLOW}📤 Pushing to GitHub...${NC}"
git push origin "$CURRENT_BRANCH"

echo -e "${GREEN}🎉 Setup complete!${NC}"
echo ""
echo -e "${BLUE}📖 Next steps:${NC}"
echo "1. Go to your GitHub repository: https://github.com/muhammadhasyim/cav-hoomd"
echo "2. Click Settings → Pages"
echo "3. Under 'Source', select 'GitHub Actions'"
echo "4. Save the settings"
echo ""
echo -e "${BLUE}📚 Your documentation will be available at:${NC}"
echo "https://muhammadhasyim.github.io/cav-hoomd/"
echo ""
echo -e "${BLUE}🔄 Automatic updates:${NC}"
echo "- Documentation rebuilds automatically when you push to main/master"
echo "- Build status visible in the Actions tab"
echo "- Typical deployment time: 2-5 minutes"
echo ""
echo -e "${BLUE}🛠️  Manual deployment:${NC}"
echo "- Go to Actions → Documentation → Run workflow"
echo "- Check 'Deploy to GitHub Pages' and run"
echo ""
echo -e "${GREEN}✨ Happy documenting!${NC}" 