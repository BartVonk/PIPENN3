from sys import argv

def load_excluded_ids(id_file_path):
    with open(id_file_path, 'r') as f:
        ids = [line.rstrip() for line in f]
        return set(ids)

def define_clusters(lines):
    clusters = []
    i = 0
    while i < len(lines) - 1:
        # Check for start of a new block
        if lines[i].startswith('>') and lines[i+1].startswith('>'):
            start = i
            i += 2
            # Continue until the next block or end of file
            while i < len(lines) - 1:
                if lines[i].startswith('>') and lines[i+1].startswith('>'):
                    break
                i += 1
            end = i if i < len(lines)-1 else len(lines)
            clusters.append((start, end))
        else:
            i += 1
    return clusters

def cluster_ids(cluster):
    unique_lines = set()
    #print('cluster:', cluster)

    for line in cluster:
        line = line.strip()
        if line.startswith('>'):
            unique_lines.add(line[1:-1])

    return set(unique_lines)

def filter_clusters(fasta_path, id_file_path, output_path):
    excluded_ids = load_excluded_ids(id_file_path)
    counter = 0

    with open(fasta_path, 'r') as infile, open(output_path, 'w') as outfile:
        input_data = list(infile)
        clusters = define_clusters(input_data)
        for cluster in clusters:
            # get ids from cluster
            ids_in_cluster = cluster_ids(input_data[cluster[0]: cluster[1]])
            #print('range:', cluster[0], cluster[1])
            #print('ids:', ids_in_cluster)
            #print()
            # if any match pass
            if excluded_ids & ids_in_cluster:
                counter+=1
                pass
            # else write cluster to output
            else:
                for line in input_data[cluster[0]: cluster[1]]:
                    outfile.write(line)
    print(counter)

if __name__ == "__main__":
    if len(argv) != 4:
        print("Usage: python script.py <input_file> <id_file> <output_file>")
        exit()
    
    input_file = argv[1]
    id_file = argv[2]
    output_file = argv[3]
    filter_clusters(input_file, id_file, output_file)
