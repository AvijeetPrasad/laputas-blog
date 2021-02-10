# Jupyter book notes

- First activate the conda environment used for the book:
```
	conda activate mybook
```
This ensures that all the packages used in the notebooks are installed in the run environment.

- If there is a major change in the notebook, run:
```
    jb clean --all book/
    jb build book/
```
Note that the "clean --all" step removes all the cache, so all the files will be executed again. Do that only if there was a major change in the file structure. For most purposes, just the "build" step would do.

- We then add and push all the files to the git hub repository:
```
    git add ./
    git commit -m "update"
    git push
    ghp-import -n -p -f book/_build/html
``` 
the last step adds these files to the github-pages branch, which is then directly published online.

For more details on the Jupyter-book visit their website: https://jupyterbook.org/intro.html

