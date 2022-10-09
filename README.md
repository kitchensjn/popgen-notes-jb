# popgen-notes-jb

/chapters - Graham's latex files
/docs - contains the static HTML from JupyterBook and used for the GitHub Pages deployment
/popgen-notes-jb - all of the source files

Use pandoc to convert latex files to markdown and then edit any errors that come up during conversion.

CMD: jupyter-book build popgen-notes-jb
creates a /_build directory within /popgen-notes-jb using the source files. I just moved /_build/html to /docs for deployment.
