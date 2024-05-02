import matplotlib.pyplot as plt
import networkx as nx
import random

class TreeNode:
    def __init__(self, name):
        self.name = name
        self.children = []
        self.maternal_cnvs = []
        self.paternal_cnvs = []
        self.maternal_fasta = None
        self.paternal_fasta = None
        self.fq1 = None
        self.fq2 = None
        self.fasta = None
        self.maternal_fasta_length = 0
        self.paternal_fasta_length = 0
        self.parent = None
        self.ratio = None
        self.cell_no = None
        self.depth = None
        self.changes = []

def draw_tree(node, pos=None, level=0, width=2., vert_gap=0.2, xcenter=0.5):
    if pos is None:
        pos = {node.name: (xcenter, 1 - level * vert_gap)}
    else:
        pos[node.name] = (xcenter, 1 - level * vert_gap)
    if node.children:
        dx = width / 2
        nextx = xcenter - width / 2 - dx / 2
        for child in node.children:
            nextx += dx
            pos = draw_tree(child, pos=pos,
                            level=level + 1,
                            width=dx, xcenter=nextx)
    return pos

def generate_random_tree(node_count):
    nodes = [TreeNode(i) for i in range(1, node_count + 1)]
    used_nodes = []
    root = random.choice(nodes)
    nodes.remove(root)
    used_nodes.append(root)
    while nodes:
        node = random.choice(nodes)
        nodes.remove(node)
        parent = random.choice(used_nodes)
        parent.children.append(node)
        node.parent = parent
        used_nodes.append(node)
    
    return root

def cal_tree_depth(tree):
    if len(tree.children) == 0:
        return 1

    return max(cal_tree_depth(child) for child in tree.children) + 1



def draw_tree_to_pdf(random_tree, outfile):
    # Creating a graph for drawing
    graph = nx.DiGraph()
    graph.add_node(random_tree.name)
    def add_edges(node):
        for child in node.children:
            graph.add_edge(node.name, child.name)
            add_edges(child)
    add_edges(random_tree)

    # Drawing the tree using matplotlib
    pos = draw_tree(random_tree)
    plt.figure(figsize=(20, 10))
    nx.draw(graph, pos, with_labels=True, arrows=False, node_size=700, node_color="skyblue", font_size=8, font_color="black")

    # Save the tree graph to a PDF file
    plt.savefig(outfile, format="pdf", dpi=300, bbox_inches="tight")

def rename_tree_nodes(root, start_value=0):
    """
    Rename tree nodes starting from the given start_value based on depth.
    
    Parameters:
    - root: The root node of the tree.
    - start_value: The starting value for renaming nodes (default is 0).
    """
    stack = [(root, 0)]
    current_value = start_value

    while stack:
        node, depth = stack.pop()

        # Rename the current node
        if depth == 0:
            node.name = 'normal'
            continue
        else:
            node.name = 'clone' + str(current_value)
        node.depth = depth
        current_value += 1

        # Add children to the stack with increased depth
        stack.extend((child, depth + 1) for child in node.children)

def generate_random_tree_balance(node_count, max_depth):
    # Example usage
    tree_depth = 1000
    while tree_depth > max_depth:
        root = generate_random_tree(node_count)
        tree_depth = cal_tree_depth(root)

    # rename node name
    rename_tree_nodes(root)
    return root

def tree_to_newick(node):
    if not node.children:
        # If the node has no children, return its name
        return node.name
    else:
        # If the node has children, recursively process each child
        child_strings = [tree_to_newick(child) for child in node.children]

        # Join child strings with commas and enclose in parentheses
        children_str = ",".join(child_strings)
        result = f"({children_str})"

        # Append the current node's name
        if node.name:
            result += node.name

        # You can add additional information if needed
        # For example, adding maternal and paternal CNVs
        # if node.maternal_cnvs:
        #     result += f"[maternal_cnvs={','.join(map(str, node.maternal_cnvs))}]"
        # if node.paternal_cnvs:
        #     result += f"[paternal_cnvs={','.join(map(str, node.paternal_cnvs))}]"
        if node.ratio:
            result += "[ratio={0}]".format(node.ratio)
        if node.cell_no:
            result += "[cell_no={0}]".format(node.cell_no)

        return result

# tree = generate_random_tree_balance(10, 4)
# draw_tree_to_pdf(tree, 'tree_graph.pdf')