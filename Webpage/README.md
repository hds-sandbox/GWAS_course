# GWAS with the sandbox workshop 

This branch hosts and generates the self-learning material for the course.

- The markdowns (*.qmd) files in `develop` and `access` are stand-alone.
- The practical lessons (`*ipynb`) are pulled from the Notebooks folder in the `main` branch. 

To keep the website up to date, stage your changes with `git add`, commit them, and push them to the remote repository. Once pushed, the `.github/workflows/publish.yml` workflow will automatically generate the updated webpage and deploy it to the `gh-pages` branch. There’s also another workflow running in the main branch that tracks changes to the notebooks (`trigger_webpage.yml`). If any of them get updated, it automatically triggers the Quarto publishing process, making sure those changes are processed and reflected on the website without any extra steps from you. If any of them get updated, it automatically triggers the Quarto publishing process, making sure those changes are processed and reflected on the website without any extra steps from you.

Check the GitHub workflows to see how this is done. 

Contributors: 
- Samuele Soraggi 
- Alba Refoyo Martinez
- Conor O'Hare
