# popgen-notes-jb

This repository hosts an online and interactive version of Graham Coop's Population And Quantitative Genetics textbook. The book has been converted using Pandoc and JupyterBook and is hosted using GitHub Pages. 

## Repository Organization

- `/chapters` - Graham's original latex files
- `/docs` - contains the static HTML from JupyterBook and used for the GitHub Pages deployment
- `/popgen-notes-jb` - all of the source files

## Updating The Book

### Set Up

`conda create --name gc-popgen-notes --file requirements.txt`

### Converting From Latex To Markdown

`pandoc -s <latex file> -o <markdown file>`

### Formatting Markdown

[Style Guide](popgen-notes-jb/formatting.md)

It may be necessary to replace `:::` with `'''` or vice versa if there are lots of nested statements. It is helpful to program these styles as shortcuts in your text editor to speed up the formatting process.

### Building HTML

`jupyter-book build popgen-notes-jb`

This creates a `/_build` directory within `/popgen-notes-jb` using the source files. You can view `/popgen-notes-jb/_build/html/index.html` in web browser to see your build.

### Deploying

Move all files in `/_build/html` to `/docs` for deployment.
