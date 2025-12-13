# routeplot.py
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

def load_xy(filename):
    """Load only first two columns, skipping lines beginning with '#'."""
    data = []
    with open(filename) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split()
            # use only first two columns as floats
            data.append([float(cols[0]), float(cols[1])])
    data.append(data[0]) # close the path from the last to the 1st city
    return np.array(data)

def load_polygons(filename):
    """Load longitude/latitude polygons separated by blank lines."""
    polygons = []
    current = []

    with open(filename) as f:
        for line in f:
            line = line.strip()

            # Blank line = end of current polygon
            if not line:
                if current:
                    polygons.append(np.array(current))
                    current = []
                continue

            # Parse longitude/latitude
            parts = line.split()
            lon = float(parts[0])
            lat = float(parts[1])
            current.append([lon, lat])

    # Add last polygon if missing trailing blank line
    if current:
        polygons.append(np.array(current))

    return polygons

def load_schedule_csv(filename):
    """
    Load annealing schedule CSV produced by tsp_sa.cpp:
      columns: T,Lbest,Lcur
    Returns: T, Lbest, Lcur (numpy arrays)
    """
    arr = np.genfromtxt(filename, delimiter=',', names=True)
    return arr['T'], arr['Lbest'], arr['Lcur']

def _infer_tag_from_paths(infile, optfile=None):
    """
    Infer tag to match required output names:
      an150.png, an1k.png, an2k.png
    We look for '150', '1k'/'1000', '2k'/'2000' in filenames.
    """
    candidates = [infile, optfile] if optfile else [infile]
    joined = " ".join([c for c in candidates if c])

    if "150" in joined:
        return "150"
    if ("1k" in joined) or ("1K" in joined) or ("1000" in joined):
        return "1k"
    if ("2k" in joined) or ("2K" in joined) or ("2000" in joined):
        return "2k"
    return None

def make_plot(infile, optfile=None, region="NA", schedule=None, save_png=False, show=True):
    '''
    infile: (required) a list of cities
    outfile: a list of cities ordered for optimized route
    region: area of the globe to plot"
            NA = North America
            World = the whole world
    schedule: optional annealing schedule CSV (schedule_TAG.csv)
    save_png: if True, also save annealing schedule as anTAG.png
    show: display the city plot window
    '''

    # --------- City route plot (original behavior) ---------
    polygons    = load_polygons("world_50m.dat")
    cities_orig = load_xy(infile)
    if optfile:
        cities_out  = load_xy(optfile)

    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot world map (outline)
    for poly in polygons:
        ax.plot(poly[:,0], poly[:,1], color="black", lw=0.8)

    # Plot original city order (thin line)
    ax.plot(cities_orig[:,0], cities_orig[:,1], lw=1, color="red", alpha=0.5)

    # Plot salesman path (lines + points)
    if optfile:
        ax.plot(cities_out[:,0], cities_out[:,1],
                lw=2, color="blue", marker='o', markersize=3)

    ax.set_title("Plot of Salesman's Cities")
    ax.set_xlabel("longitude")
    ax.set_ylabel("latitude")

    # Axis ranges
    ax.set_xlim(-180, -60)
    ax.set_ylim(10, 75)
    if region == "World":
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Save route plot (unchanged default: pdf)
    plotfile = infile.split('.')[0] + '.pdf'
    plt.savefig(plotfile, format='pdf', facecolor='white')

    if schedule:
        T, Lbest, Lcur = load_schedule_csv(schedule)

        fig2, ax2 = plt.subplots(figsize=(7.5, 5.0))
        ax2.plot(T, Lbest, label="best")
        ax2.plot(T, Lcur,  label="current", alpha=0.6)

        ax2.set_xscale("log")
        ax2.set_yscale('log')
        ax2.invert_xaxis()  # cooling: high T -> low T

        ax2.set_xlabel("Temperature T")
        ax2.set_ylabel("Total distance (km)")
        ax2.set_title("Annealing Schedule: Distance vs Temperature")
        ax2.legend()

        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)

        # Save schedule plot: required names an150.png / an1k.png / an2k.png
        if save_png:
            tag = _infer_tag_from_paths(infile, optfile)
            if tag is None:
                # fallback: derive from schedule filename
                base = os.path.basename(schedule)
                # schedule_150.csv -> 150, schedule_1k.csv -> 1k, schedule_2k.csv -> 2k
                tag = base.replace("schedule_", "").replace(".csv", "")
            outpng = f"an{tag}.png"
            fig2.savefig(outpng, dpi=200, facecolor='white')
            print(f"Saved annealing schedule plot: {outpng}")

    if show:
        print('close plot or "^C" to exit')
        try:
            plt.show()
        except KeyboardInterrupt:
            print("Interrupted with Ctrl-C, closing plot and exiting...")
            plt.close('all')


        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='route plotter')
    parser.add_argument('paths', nargs='*',
                        help='''paths to plot: original [optimized]''')
    parser.add_argument("-w", action='store_true',
                        help="plot the whole world, default in North America only")
    # NEW args
    parser.add_argument("-s", "--schedule", default=None,
                        help="annealing schedule CSV (e.g., schedule_150.csv)")
    parser.add_argument("--save-png", action='store_true',
                        help="save annealing schedule plot as anTAG.png (an150.png/an1k.png/an2k.png)")
    parser.add_argument("--no-show", action='store_true',
                        help="do not open interactive windows (just save files)")

    args = parser.parse_args()

    if len(args.paths) < 1:
        print("at least one input file needed")
        exit(1)

    cities = args.paths[0]
    cities2 = args.paths[1] if len(args.paths) > 1 else None
    region = "World" if args.w else "NA"

    make_plot(
        cities,
        cities2,
        region,
        schedule=args.schedule,
        save_png=args.save_png,
        show=(not args.no_show),
    )
