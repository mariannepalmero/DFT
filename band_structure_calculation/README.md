# Band Structure Calculation Using Quantum ESPRESSO

Here are the steps required to calculate a simple band structure of an optimized system using Quantum ESPRESSO.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Step 1: Conduct the Self-Consistent Field (SCF) Calculation
> Perform an SCF calculation on your optimized system.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Step 2: Prepare for Band Structure Calculation
> 1. Create a new folder named "bands" in your working directory.
> 2. Copy the INPUT file from the SCF calculation into the bands folder.
> 3. Modify the following lines in the copied INPUT file: \
> a. Change the calculation type from `scf` to `bands`. \
> b. Update the outdir path from `./tmp/` to `../tmp/`, to ensure the correct output directory is used. \
> c. Modify the `K_POINTS {angstrom}` section to `K_POINTS crystal` format. 

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Step 3: Obtain the K_POINTS crystal Lines
> To generate the correct K_POINTS crystal lines, follow these steps:
> 
> 1. Open the output file from the SCF calculation using XCrySDen.
> 2. Click on Tools → K-Path Selection.
> 3. In the Primitive Brillouin Zone window, select the main path by clicking on the dots for the high-symmetry points (e.g., Γ → M → K → Γ).
> 4. Adjust the number of k-points along the path to your desired resolution (e.g., 100 points).
> 5. Save the file with the .pwscf extension.
> 6. Open the saved file and copy the text corresponding to the K_POINTS crystal lines into your input file.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Step 4: Run the Band Structure Calculation
> Execute the band structure calculation using the pw.x command. Ensure that the proper input file with the bands calculation setting is specified (eg., `pw.x -in bands_input.in > bands_output.out`)

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Step 5: Post-Processing the Band Structure Data
> 1. Create a new folder named "pp" for post-processing.
> 2. Generate an INPUT file for the bands.x post-processing utility. A sample input file might look like this:
<pre>
&BANDS
 prefix = 'surf-1' 
 outdir = '../../tmp' 
 filband = 'surf-1.dat' 
/
</pre>
> If you are calculating spin-polarized systems, create separate pp folders for spin-up and spin-down calculations. Add the line `spin_component = 1` (for spin-up) or `spin_component = 2` (for spin-down) in the respective input files.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Step 6: Run the Post-Processing Calculation
> Run the bands.x utility with the prepared input file (eg., `bands.x -in bands_pp.in > bands_pp.out`).

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Step 7: Plot the Band Structure
> Once the band structure data is ready, use the `plot_bands.py` to visualize the results. You can modify the script to suit your plotting preferences, such as adjusting labels, colors, or axis limits (eg., `python plot_bands.py`).

<p align="right">(<a href="#readme-top">back to top</a>)</p>
