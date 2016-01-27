import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Scanner;
import java.util.TreeMap;

public class PhylogeneticTree {
	
	public static class Node {
		public String name, parName;
		public double time = 0;
		public List<Node> children = new ArrayList<Node>(2);
		public int index = -1;
		public boolean isLeaf(){ return children.size()==0; }
	}
	
	/** Array of nodes used for easy traversing of tree (leaves->root) */
	private Node[] arr;
	
	/**
	 * Loads a phylogenetic tree from file. Input file has lines in form:
	 * "[node name] [parent name] [time]". Lines can be in any order (they
	 * don`t have to follow root->leaves or leaves->root order).
	 */
	public PhylogeneticTree(String filePath){
		// Read input file, convert to map [name]->[node object]
		Map<String, Node> map = new TreeMap<String, Node>();
		Node root = new Node();
		root.name = "Root";
		map.put("Root", root);
		try {
			Scanner in = new Scanner(new File(filePath));
			in.useLocale(Locale.US);
			while(in.hasNext()){
				Node n = new Node();
				n.name = in.next();
				n.parName = in.next();
				n.time = in.nextDouble();
				map.put(n.name, n);
			}
			// Link parent-child nodes
			for(Map.Entry<String, Node> e:map.entrySet()){
				Node n = e.getValue();
				if(n == root){ continue; }
				Node par = map.get(n.parName);
				par.children.add(n);
			}
			in.close();
		}
		catch (FileNotFoundException e) {e.printStackTrace();}
		
		// Convert map to array using breadth-first search
		arr = new Node[map.size()];
		int i = map.size()-1, freeI = i-1;
		arr[i] = root;
		for(; i >= 0; i--){
			Node n = arr[i];
			n.index = i;	// The node object needs to know its index
			for(Node c:n.children){
				arr[freeI--] = c;
			}
		}
	}
	
	/** Returns node of specified array index */
	public Node get(int i){
		if(i<0 || i>arr.length-1){ return null; }
		return arr[i];
	}
	
	/** Returns total number of nodes (not only leaves) */
	public int size(){
		return arr.length;
	}
}