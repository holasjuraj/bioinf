import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;

public class Felsenstein {
	public static final char[] bases = {'A', 'C', 'T', 'G'};
	
	/**
	 * Implementation of Felsenstein algorithm to compute probability of given
	 * phylogenetic tree with data in its leaves. Method uses and Jukes-Cantor
	 * evolution model.
	 * @param tree : phylogenetic tree
	 * @param alpha : evolution speed
	 * @param data : data in leaves (map [animal]->[base])
	 */
	public static double treeProbability(PhylogeneticTree tree, double alpha, Map<String, Character> data){
		double[][] A = new double[tree.size()][4];
		
		for(int nodeI = 0; nodeI < tree.size(); nodeI++){
			PhylogeneticTree.Node n = tree.get(nodeI);
			for(int baseI = 0; baseI < 4; baseI++){
				if(n.isLeaf()){
					// Leaf
					char dataBase = data.get(n.name).charValue();
					if(dataBase==bases[baseI] || dataBase=='N' || dataBase=='-'){
						A[nodeI][baseI] = 1;
					}
					else{
						A[nodeI][baseI] = 0;						
					}
				}
				else{
					// Inner node
					PhylogeneticTree.Node childL = n.children.get(0),
										  childR = n.children.get(1);
					double sumL = 0, sumR = 0;
					for(int chBaseI = 0; chBaseI < 4; chBaseI++){
						sumL += modelJC(baseI, chBaseI, childL.time, alpha) * A[childL.index][chBaseI];
						sumR += modelJC(baseI, chBaseI, childR.time, alpha) * A[childR.index][chBaseI];
					}
					A[nodeI][baseI] = sumL * sumR;
				}
			}
		}
		
		// Root
		double sum = 0;
		for(int baseI = 0; baseI < 4; baseI++){
			sum += A[tree.size()-1][baseI] * 0.25;	// equilibrium probabilities q_a = 0.25
		}
		
		return sum;
	}
	
	/**
	 * Finds optimal value of alpha for given sequence alignment. Method uses
	 * Felsenstein algorithm and Jukes-Cantor evolution model. Alpha is from
	 * range 0:0.1:2.
	 * @param seqs : aligned sequences (map [animal]->[sequence]). All
	 * 				 sequences must be same length.
	 * @param tree : phylogenetic tree
	 */
	public static double bestAlpha(TreeMap<String, String> seqs, PhylogeneticTree tree){
		int seqLength = seqs.get("Human").length();
		double bestAlpha = -1;
		LogNum bestProb = new LogNum(0);
		Map<String, Character> column = new HashMap<>();
		
		for(double alpha = 0; alpha <= 2.0001; alpha += 0.1){
			// Compute alpha probability
			LogNum alphaProb = new LogNum(1);
			for(int i = 0; i < seqLength; i++){
				// Create column data
				for(Iterator<Entry<String, String>> it = seqs.entrySet().iterator(); it.hasNext();){
					Entry<String, String> e = it.next();
					column.put(e.getKey(), e.getValue().charAt(i));
				}
				// Compute column probability
				double prob = treeProbability(tree, alpha, column);
				alphaProb = LogNum.mul(prob, alphaProb);
			}
			
			// Compare to best
			if(alphaProb.isGt(bestProb)){
				bestProb = alphaProb;
				bestAlpha = alpha;
			}
		}
		
		return bestAlpha;
	}
	
	/**
	 * Computes probability of changing oldBase to newBase in time t using
	 * Jukes-Cantor evolution model.
	 */
	public static double modelJC(int oldBase, int newBase, double t, double alpha){
		// Jukes-Cantor evolution model
		if(oldBase == newBase){
			return (1 + 3*Math.exp(-4/3*alpha*t))/4;
		}
		else{
			return (1 - Math.exp(-4/3*alpha*t))/4;
		}
	}
	
}