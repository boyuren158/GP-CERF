# Reproducing Software Paper Examples

Follow these steps to reproduce the software paper examples:

* Step 1: Clone the repo and navigate to the folder.

```
git clone -b JOSS https://github.com/NSAPH-Software/GPCERF.git
cd GPCERF
```


* Step 2: Run docker container. In the terminal, run the following:

```
docker run -it --rm \
        -p 8230:8787 \
        -e USER=rstudio \
        -e PASSWORD=pass \
        -v "/path/to/your/folder/on/host:/home/rstudio/Project" nsaphsoftware/gpcerf_dev:3
```

* Step 3: Open your browser and go to the following link:

```
http://localhost:8230/
```

* Step 4: Log in with the following credentials:

- Username: `rstudio`
- Password: `pass`

* Step 5: Load the project and run the examples:

```
devtools::load_all()
devtools::install()
# and run the examples.
```
