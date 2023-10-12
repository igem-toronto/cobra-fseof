<div align="center">    

# cobra-fseof

</div>

<p align="center">
 <img width="450" src="https://www.ncbi.nlm.nih.gov/corecgi/tileshop/tileshop.fcgi?p=PMC3&id=449219&s=21&r=1&c=1">
</p>

cobra-fseof is a COBRApy implementation of the FSEOF (**F**lux **S**canning based on **E**nforced **O**bjective **F**lux) 
algorithm [(Choi et al., 2010)](https://doi.org/10.1128/AEM.00115-10). The image above is taken from Figure 1 of their paper. 
This library supports one core function: 

```python
def fseof(
    model: cobra.Model,
    num_steps: int,
    enforced: ReactionLike,
    objective: ReactionLike,
    enforced_direction: Literal["min", "max"] = "max",
    enforced_frac_opt: float = 0.9,
    reactions: Optional[List[ReactionLike]] = None,
    optimizer: Optional[Callable[[cobra.Model], cobra.Solution]] = None,
) -> FSEOFResults:
```

which can be used as follows:

```python
import cobra_fseof

results = cobra_fseof.fseof(...)
```

Please refer to the source code for further documentation.

______________________________________________________________________

## Installation

Pip from source:

```bash
git clone https://github.com/igem-toronto/cobra-fseof
 
cd cobra-fseof
pip install -e .   
 ```

______________________________________________________________________

## Usage 

### Lycopene Production 
 
In the `examples/lycopene` folder, we try to reproduce the results from Choi et al. (2010),
where they use FSEOF to find gene amplification targets to enhance lycopene production in *E. coli*.

### iGEM Toronto 2023

This repository was created as a part of iGEM Toronto 2023 project.
For further details, we refer readers to the iGEM [homepage](https://igem.org/)
and our team [wiki](https://2023.igem.wiki/toronto/). 