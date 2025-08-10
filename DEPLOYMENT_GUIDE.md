# 📚 GitHub Pages Deployment Guide for splikit

## Overview
This guide will help you deploy the pkgdown documentation site to GitHub Pages.

## 🚀 Deployment Steps

### Step 1: Push Changes to GitHub
```bash
# Add all changes
git add .

# Commit (if you haven't already)
git commit -m "Ready for pkgdown deployment"

# Push to your branch
git push origin version_2

# If version_2 is not your default branch, merge to main
git checkout main
git merge version_2
git push origin main
```

### Step 2: Enable GitHub Actions (if not already enabled)
1. Go to your repository: https://github.com/Arshammik/splikit
2. Click on **Actions** tab
3. If prompted, enable GitHub Actions for the repository

### Step 3: Configure GitHub Pages
1. Go to **Settings** → **Pages** (left sidebar)
2. Under **Source**, select: **Deploy from a branch**
3. Under **Branch**, select: **gh-pages**
4. Under **Folder**, select: **/ (root)**
5. Click **Save**

### Step 4: Trigger the Build
The pkgdown site will build automatically when you:
- Push to main/master branch
- Create a release
- Manually trigger (Actions tab → pkgdown → Run workflow)

To manually trigger:
1. Go to **Actions** tab
2. Select **pkgdown** workflow
3. Click **Run workflow** → **Run workflow**

### Step 5: Monitor the Build
1. Go to **Actions** tab
2. Click on the running workflow
3. Watch the progress (takes 3-5 minutes)
4. ✅ Green checkmark = Success!

### Step 6: Access Your Site
Once deployed, your site will be available at:
**https://arshammik.github.io/splikit/**

First deployment may take 5-10 minutes to become visible.

## 🔧 Troubleshooting

### If the build fails:
1. Check the Actions tab for error messages
2. Common issues:
   - **"No such file or directory"**: The .nojekyll file should fix this
   - **"Articles not found"**: We've removed articles section
   - **"Package won't install"**: Check DESCRIPTION dependencies

### If the site doesn't appear:
1. Check Settings → Pages to ensure it's enabled
2. Wait 10 minutes (first deployment can be slow)
3. Try: https://arshammik.github.io/splikit/index.html
4. Clear browser cache (Ctrl+Shift+R or Cmd+Shift+R)

### To rebuild the site:
```bash
# Make changes to documentation
# Then commit and push
git add .
git commit -m "Update documentation"
git push origin main

# The GitHub Action will automatically rebuild
```

## 📁 File Structure

Your repository now has:
```
splikit/
├── .github/
│   └── workflows/
│       └── pkgdown.yaml    # Builds and deploys site
├── _pkgdown.yml            # Site configuration
├── .nojekyll              # Prevents Jekyll processing
├── docs/                  # Will contain built site (on gh-pages branch)
└── vignettes/             # Your documentation
```

## 🎨 Customization

To customize the site appearance, edit `_pkgdown.yml`:
```yaml
template:
  bootstrap: 5
  bootswatch: cosmo  # Try: cerulean, cosmo, flatly, journal, etc.
  theme: github-light  # or github-dark
```

## 📋 Checklist

- [ ] Push all changes to GitHub
- [ ] GitHub Actions is enabled
- [ ] pkgdown workflow runs successfully
- [ ] GitHub Pages is configured for gh-pages branch
- [ ] Site is accessible at https://arshammik.github.io/splikit/

## 🔗 Useful Links

- **Your Site**: https://arshammik.github.io/splikit/
- **Repository**: https://github.com/Arshammik/splikit
- **Actions**: https://github.com/Arshammik/splikit/actions
- **Settings**: https://github.com/Arshammik/splikit/settings/pages

## 💡 Tips

1. **Local Preview**: Build locally before pushing
   ```r
   pkgdown::build_site()
   # Open docs/index.html in browser
   ```

2. **Update Only Reference**: For quick updates
   ```r
   pkgdown::build_reference()
   ```

3. **Check Configuration**:
   ```r
   pkgdown::check_pkgdown()
   ```

## 📝 Notes

- The site rebuilds automatically on every push to main
- The gh-pages branch is managed automatically by GitHub Actions
- Don't edit the gh-pages branch directly
- The .nojekyll file ensures proper deployment without Jekyll

---

Good luck with your deployment! 🎉