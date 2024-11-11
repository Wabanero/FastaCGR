import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import argparse
import os
import logging
from multiprocessing import Pool, cpu_count
from functools import partial

def create_directories(directories):
    """
    Creates the necessary directories if they don't exist.
    """
    for directory in directories:
        try:
            os.makedirs(directory, exist_ok=True)
            print(f"Directory '{directory}' is ready.")
        except Exception as e:
            print(f"Error creating directory '{directory}': {e}")

def configure_logging(log_dir, log_file):
    """
    Configures the log settings.
    """
    log_path = os.path.join(log_dir, log_file)
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logging.info("Logging configured successfully.")

def log_message(message):
    """
    Logs a message to both console and log file.
    """
    print(message)
    logging.info(message)

def read_fasta(file_path):
    """
    Reads a FASTA file and extracts sequences.
    """
    sequences = []
    try:
        with open(file_path, 'r') as f:
            header = None
            seq = ''
            for line in f:
                if line.startswith('>'):
                    if header and seq:
                        sequences.append((header, seq))
                        seq = ''
                    header = line[1:].strip()
                else:
                    seq += line.strip()
            if header and seq:
                sequences.append((header, seq))
    except FileNotFoundError:
        log_message(f"Error: File {file_path} not found.")
    except Exception as e:
        log_message(f"Error reading FASTA file: {e}")
    return sequences

def generate_cgr(sequence, dimension=2, k=1):
    """
    Generates a CGR or FCGR matrix for a given sequence.
    If k > 1, it generates an FCGR by counting k-mer frequencies.

    Parameters:
        sequence (str): The nucleotide sequence.
        dimension (int): Dimension of the CGR (2 or 3).
        k (int): Length of k-mers for FCGR. If k=1, standard CGR is generated.

    Returns:
        numpy.ndarray: CGR or FCGR matrix.
    """
    if dimension not in [2, 3]:
        raise ValueError("Dimension must be 2 or 3.")

    sequence = sequence.upper()

    # Define the vertices for nucleotides
    if dimension == 2:
        vertices = {'A': [0, 0], 'C': [0, 1], 'G': [1, 1], 'T': [1, 0]}
    else:
        # Map nucleotides to vertices of a tetrahedron
        vertices = {
            'A': [1, 1, 1],
            'C': [1, 0, 0],
            'G': [0, 1, 0],
            'T': [0, 0, 1],
            'U': [0, 0, 1]  # U is mapped to T for RNA sequences
        }

    # Initialize the starting point at the center
    if dimension == 2:
        point = [0.5, 0.5]
    else:
        point = [0.5, 0.5, 0.5]
    cgr_points = [point.copy()]

    if k == 1:
        # Generate standard CGR
        for nucleotide in sequence:
            if nucleotide in vertices:
                vertex = vertices[nucleotide]
                point = [(p + v) / 2 for p, v in zip(point, vertex)]
                cgr_points.append(point.copy())
            else:
                # Skip invalid characters
                continue
        cgr_matrix = np.array(cgr_points)
    else:
        # Generate FCGR
        num_bins = 2 ** k
        if dimension == 2:
            fcgr_matrix = np.zeros((num_bins, num_bins))
        else:
            fcgr_matrix = np.zeros((num_bins, num_bins, num_bins))

        total_kmers = len(sequence) - k + 1
        for i in range(total_kmers):
            kmer = sequence[i:i+k]
            if all(nuc in vertices for nuc in kmer):
                position = point.copy()
                for nuc in kmer:
                    vertex = vertices[nuc]
                    position = [(p + v) / 2 for p, v in zip(position, vertex)]
                # Map position to grid index
                indices = (np.array(position) * (num_bins - 1)).astype(int)
                indices = tuple(indices)
                fcgr_matrix[indices] += 1
        cgr_matrix = fcgr_matrix

    return cgr_matrix

def plot_cgr(cgr_matrix, output_file, dimension=2, cmap='viridis', width=800, height=800):
    """
    Plots and saves the CGR or FCGR image.
    """
    if dimension == 2:
        plt.figure(figsize=(width/100, height/100))
        if cgr_matrix.ndim == 2 and cgr_matrix.shape[1] == 2:
            # Standard CGR plot
            plt.scatter(cgr_matrix[:, 0], cgr_matrix[:, 1], s=0.1, c='black')
        else:
            # FCGR heatmap
            plt.imshow(cgr_matrix, cmap=cmap, origin='lower')
        plt.axis('off')
        plt.savefig(output_file, bbox_inches='tight', pad_inches=0, dpi=100)
        plt.close()
    else:
        fig = plt.figure(figsize=(width/100, height/100))
        ax = fig.add_subplot(111, projection='3d')
        if cgr_matrix.ndim == 2 and cgr_matrix.shape[1] == 3:
            # Standard 3D CGR plot
            ax.scatter(cgr_matrix[:, 0], cgr_matrix[:, 1], cgr_matrix[:, 2],
                       s=0.5, c='black')
        else:
            # FCGR 3D heatmap
            indices = np.argwhere(cgr_matrix > 0)
            counts = cgr_matrix[cgr_matrix > 0]
            ax.scatter(indices[:, 0], indices[:, 1], indices[:, 2],
                       s=counts, c=counts, cmap=cmap)
        # Set axis labels
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        # Optionally, adjust axis properties
        ax.grid(True)
        ax.set_xticks([0, 0.5, 1])
        ax.set_yticks([0, 0.5, 1])
        ax.set_zticks([0, 0.5, 1])
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])
        ax.set_zlim([0, 1])
        plt.savefig(output_file, bbox_inches='tight', pad_inches=0, dpi=100)
        plt.close()

def process_sequence(idx, header, sequence, args):
    """
    Processes a single sequence to generate and save CGR or FCGR images.
    """
    try:
        log_message(f"[Sequence {idx}] Header: {header}")
        log_message(f"[Sequence {idx}] Length: {len(sequence)}")

        # Determine the file prefix based on k-mer length
        if args.kmer == 1:
            prefix = 'cgr'
        else:
            prefix = f'fcgr_k{args.kmer}'

        # Generate 2D CGR or FCGR
        cgr_2d = generate_cgr(sequence, dimension=2, k=args.kmer)
        output_file_2d = os.path.join(args.output, f'{prefix}_2d_{idx}.png')
        plot_cgr(cgr_2d, output_file_2d, dimension=2, cmap=args.cmap,
                 width=args.width, height=args.height)
        log_message(f"[Sequence {idx}] Saved 2D image to {output_file_2d}")

        # Generate 3D CGR or FCGR
        cgr_3d = generate_cgr(sequence, dimension=3, k=args.kmer)
        output_file_3d = os.path.join(args.output, f'{prefix}_3d_{idx}.png')
        plot_cgr(cgr_3d, output_file_3d, dimension=3, cmap=args.cmap,
                 width=args.width, height=args.height)
        log_message(f"[Sequence {idx}] Saved 3D image to {output_file_3d}")

    except Exception as e:
        log_message(f"[Sequence {idx}] Error processing sequence: {e}")

def main():
    """
    Main function.
    """
    directories = ["logs", "plots", "results", "output/images"]
    create_directories(directories)

    configure_logging("logs", "cgr_generator.log")

    # Parameters
    parser = argparse.ArgumentParser(
        description='Generate 2D and 3D CGR or FCGR plots from a multi-FASTA file.'
    )
    parser.add_argument('input', help='Input FASTA file')
    parser.add_argument('-o', '--output', help='Output directory', default='output/images')
    parser.add_argument('--width', type=int, help='Image width in pixels', default=800)
    parser.add_argument('--height', type=int, help='Image height in pixels', default=800)
    parser.add_argument('--kmer', type=int, help='Length of k-mers for FCGR (default: 1 for CGR)', default=1)
    parser.add_argument('--cmap', type=str, help='Colormap for CGR plots (default: viridis)', default='viridis')
    parser.add_argument('--processes', type=int, help='Number of parallel processes (default: number of CPU cores)', default=cpu_count())
    args = parser.parse_args()

    # Logging
    log_message(f"Reading sequences from {args.input}")
    sequences = read_fasta(args.input)
    if not sequences:
        log_message("No sequences found. Exiting.")
        return

    if not os.path.exists(args.output):
        try:
            os.makedirs(args.output)
            log_message(f"Created output directory at {args.output}")
        except Exception as e:
            log_message(f"Error creating output directory: {e}")
            return

    log_message(f"Starting processing with {args.processes} parallel processes.")

    pool = Pool(processes=args.processes)
    process_func = partial(process_sequence, args=args)
    tasks = [(idx + 1, header, seq) for idx, (header, seq) in enumerate(sequences)]
    try:
        pool.starmap(process_func, tasks)
    except KeyboardInterrupt:
        log_message("Processing interrupted by user.")
    except Exception as e:
        log_message(f"An error occurred during multiprocessing: {e}")
    finally:
        pool.close()
        pool.join()
        log_message("Processing completed.")

if __name__ == '__main__':
    main()
