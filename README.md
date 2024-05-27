
<h1 align="center">
MOLInterface
</h1>

<br>


Repository for charing codes regarding the semester project

## 👨‍🔬 Usage

> This is a program that will take as an input the current clipboard of the user, as long as it is a molecule's name.
> It will give a lot of useful infos about the molecule, like its IUPAC name, SMILES code, atomic formula, the molar mass, and will draw the molecule in 2D in a tkinter window.
> If stereogenic centers are present in the molecule, they will be given with their position and configuration (R,S) and will be displayed on the molecule.
> The chemical groups present in the molecule will also be given, as well as the rings found with their aromaticity.
> A button is also available in order to highlight the chemical groups of the molecule
> A link leading to the PUBchem page of the selected molecule is also available at the bottom of the window.

## 👇👾 Installation

Create a new environment, you may also give the environment a different name. 

```
conda create -n MOLInterface python=3.10 
```

```
conda activate MOLInterface
```

If you need jupyter lab, install it 

```
(MOLInterface) $ pip install jupyterlab
```
Finally, create a clone of the repository of your local machine in order to have access to the python package

```
git clone https://github.com/EAJ7/MOLInterface
```

## 🛠️ Development installation

Initialize Git (only for the first time). 

Note: You should have create an empty repository on your github.

```
git init 
git add .
git commit -m "Initial commit" 
git branch -M main
git remote add origin git@github.com:(your_github_profile)/(YOUR_REPOSITORY).git 
git push -u origin main
```

Then add and commit changes as usual. 

To install the package, run

```
(MOLInterface) $ pip install -e ".[test,doc]"
```
## Code usage and accessibilities

After successfuly cloning the repository into your local machine, you should have access to the code via the jupyter notebook MOLInterface.ipynb or the python file MOLInterface.py.
In addition, your newly created conda environment has access to the python package molinterface, which contains all of the functions used for running the code.
You can import these functions into your python code by applying the following code lines:

```
from molinterface import THE_FUNCTION_THAT_YOU_WANT_TO_IMPORT
```
Alternatively, if you plan on using multiple functions contained in the package, you can directly import it into your code:

```
import molinterface as moli
```
Note that some functions of this code will take the smarts_patterns dictionnary as an argument. This dictionnary is defined as a global variable in the main() function.

### Run tests and coverage
In ordrer to check that your code works properly, you can run the following commands on your terminal.
```
(MOLInterface) $ pip install tox
(MOLInterface) $ tox
```

### Generate coverage badge

Works after running `tox`

```
(MOLInterface) $ pip install "genbadge[coverage]"
(MOLInterface) $ genbadge coverage -i coverage.xml
```

Generated with some inspiration from [cookiecutter-snekpack](https://github.com/cthoyt/cookiecutter-snekpack) and [copier-pylib](https://github.com/astrojuanlu/copier-pylib).

