import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;

public class Main {

	public static class ExonPos{
		public int start, end;
		public ExonPos(int start, int end){
			this.start = start;
			this.end = end;
		}
	}
	
	public static void main(String[] args) {
		PhylogeneticTree tree = new PhylogeneticTree("data/strom.txt");
		String[] animals = {"Human", "Chimp", "Baboon", "Macaque", "Marmoset",
							"Rat", "Mouse", "Rabbit", "Horse", "Cat", "Dog",
							"Cow", "Pig", "Elephant"};
		
		/**
		 * b)
		 * Create dummy data - same base for all animals in tree, and compute
		 * tree probability with Felsenstein algorithm for 2 different values
		 * of alpha.
		 */
		Map<String, Character> data = new HashMap<>();
		for(String animal:animals){
			data.put(animal, 'A');
		}
		System.out.println("b)\nalpha = 10^9 -> Pr=" + Felsenstein.treeProbability(tree, 1000000000, data));
		System.out.println("alpha = 1.0 -> Pr=" + Felsenstein.treeProbability(tree, 1, data));
		System.out.println("alpha = 0.2 -> Pr=" + Felsenstein.treeProbability(tree, 0.2, data));
		System.out.println("alpha = 0.0 -> Pr=" + Felsenstein.treeProbability(tree, 0, data));
		
		
		/**
		 * d)
		 * Split aligned sequences into windows of length w and find best
		 * value of alpha for each window. Print result for first 10 windows.
		 */
		TreeMap<String, String> wholeSeqs = readSequences("data/cftr.txt");
		List<Double> bestAlphas = new ArrayList<>();
		int w = 100,
			seqLength = wholeSeqs.get("Human").length();
		
		for(int i = 0; i*w < seqLength; i++){
			// Prepare window
			TreeMap<String, String> window = new TreeMap<>();
			for(Iterator<Entry<String, String>> it = wholeSeqs.entrySet().iterator(); it.hasNext();){
				Entry<String, String> e = it.next();
				int start = i*w, end = Math.min((i+1)*w, seqLength);
				window.put(e.getKey(), e.getValue().substring(start, end));
			}
			// Find best alpha
			bestAlphas.add(Felsenstein.bestAlpha(window, tree));
		}
		
		// Print first 10
		System.out.println("d)\nwindow\tbest alpha");
		for(int i = 0; i < 10; i++){
			System.out.println((i+1) + "\t" + bestAlphas.get(i));
		}
		
		
		/**
		 * e)
		 * Find windows which overlap with exons in Human genome. Compute
		 * histogram for windows with exons, windows without exons, and all
		 * windows.
		 */
		String human = wholeSeqs.get("Human");
		
		// Read and save exon positions from file, change indexing to 0-bsaed
		List<Integer> exons = new LinkedList<>();		
		try {
			Scanner in = new Scanner(new File("data/exony.txt"));
			while(in.hasNext()){
				exons.add(in.nextInt()-1);	// start of exon
				exons.add(in.nextInt()-1);	// end of exon
			}
			exons.add(Integer.MAX_VALUE);	// dummy last exon
			exons.add(Integer.MAX_VALUE);
			in.close();
		}
		catch (FileNotFoundException e) { e.printStackTrace(); }
		
		// Find windows with exons
		Set<Integer> withExons = new HashSet<>();
		Iterator<Integer> exIt = exons.iterator();
		int nextStart = exIt.next(), nextEnd = exIt.next();
		boolean inExon = false, exonInWin = false;
		
		for(int i = 0, j = 0; i < seqLength; i++){
			// i is 0-based index in aligned sequence
			// j is 0-based index in dash-free sequence
			if(human.charAt(i) != '-'){
				if(!inExon && nextStart==j){	// start of exon
					inExon = exonInWin = true;
				}
				else if(inExon && nextEnd==j){		// end of exon
					inExon = false;
					nextStart = exIt.next();
					nextEnd = exIt.next();
				}
				j++;
			}
			
			// Window boundary
			if(i%w == w-1 || i == seqLength-1){
				if(exonInWin){
					withExons.add(i/w);
					// If exon ended, reset exonInWin for next window,
					// if exon overlaps to next window, keep exonInWin true 
					exonInWin = inExon;
				}
			}
		}
			
		// Count windows with/without exons for different values of alpha
		int count[][] = new int[2][21];
		for(int i = 0; i < bestAlphas.size(); i++){
			if(withExons.contains(i)){
				count[0][(int) Math.round(bestAlphas.get(i) * 10)]++;
			}
			else{
				count[1][(int) Math.round(bestAlphas.get(i) * 10)]++;
			}
		}
		double  exSize = withExons.size(),
				noExSize = bestAlphas.size() - exSize;
		
		// Compute and print frequencies
		System.out.println( "e)\n"+
							"windows with exons: " + exSize +"\n"+
							"windows without exons: " + noExSize +"\n"+
							"total windows: " + (exSize+noExSize) +"\n"+
							"frequencies:\n"+
							"with\twithout\ttotal");
		for(int a = 0; a <= 20; a++){
			System.out.println( ((double)count[0][a] / exSize) + "\t" +
								((double)count[1][a] / noExSize) + "\t" +
								((double)(count[0][a]+count[1][a]) / (exSize+noExSize)) );
		}
		
		
		/**
		 * f)
		 * Using Hidden Markov Model to locate exons in human genome. First 70%
		 * of genome is used for training, remaining 30% for testing. Genome is
		 * processed by non-overlapping windows of fixed size. Computed alphas
		 * for windows are saved to file and loaded again for quicker testing
		 * of HMM (just comment out the section with computing alphas).
		 */
		List<Integer> bestAlphasHMM = new ArrayList<>();
		int winSize = 10;
		seqLength = 0;
		
		// Filter only bases present in human genome
		TreeMap<String, StringBuilder> builders = new TreeMap<>();
		for(Iterator<Entry<String, String>> it = wholeSeqs.entrySet().iterator(); it.hasNext();){
			builders.put(it.next().getKey(), new StringBuilder());
		}
		for(int i = 0; i < human.length(); i++){
			if(human.charAt(i) == '-'){ continue; }
			for(Iterator<Entry<String, StringBuilder>> it = builders.entrySet().iterator(); it.hasNext();){
				Entry<String, StringBuilder> e = it.next();
				e.getValue().append(wholeSeqs.get(e.getKey()).charAt(i));
			}
			seqLength++;
		}
		TreeMap<String, String> filteredSeqs = new TreeMap<>();
		for(Iterator<Entry<String, StringBuilder>> it = builders.entrySet().iterator(); it.hasNext();){
			Entry<String, StringBuilder> e = it.next();
			filteredSeqs.put(e.getKey(), e.getValue().toString());
		}
		
		// Compute alphas and write to file
		try {
			PrintWriter pw = new PrintWriter(new File("data/alphas.txt"));
			
			for(int i = 0; i*winSize < seqLength; i++){
				// Prepare window
				TreeMap<String, String> window = new TreeMap<>();
				for(Iterator<Entry<String, String>> it = filteredSeqs.entrySet().iterator(); it.hasNext();){
					Entry<String, String> e = it.next();
					int start = i*winSize, end = Math.min((i+1)*winSize, seqLength);
					window.put(e.getKey(), e.getValue().substring(start, end));
				}
				// Find best alpha
				pw.println((int) Math.round(Felsenstein.bestAlpha(window, tree) * 10));
			}
			
			pw.flush();
			pw.close();
		}
		catch (IOException e) { e.printStackTrace(); }
		
		// Read individual alphas from file 
		try {
			Scanner in = new Scanner(new File("data/alphas.txt"));
			while(in.hasNext()){
				bestAlphasHMM.add(in.nextInt());
			}
			in.close();
		}
		catch (FileNotFoundException e) { e.printStackTrace(); }
		
		// Prepare label string (HMM states) - find windows with exons
		StringBuilder sb = new StringBuilder();
		Iterator<Integer> it = exons.iterator();
		nextStart = it.next(); nextEnd = it.next();
		for(int i = 0; i < bestAlphasHMM.size(); i++){
			if(nextStart >= (i+1)*winSize){
				sb.append(HMM.state1);
			}
			else{
				sb.append(HMM.state2);
				while(nextEnd < (i+1)*winSize){
					nextStart = it.next(); nextEnd = it.next();
				}
			}
		}
		String label = sb.toString();
		
		// Training and testing HMM
		int trainSize = (bestAlphasHMM.size() * 70) / 100;
		List<Integer> trainAlphas = bestAlphasHMM.subList(0, trainSize);
		String trainLabel = label.substring(0, trainSize);
		List<Integer> testAlphas = bestAlphasHMM.subList(trainSize, bestAlphasHMM.size());
		String testLabel = label.substring(trainSize, label.length());
		
		System.out.println("f) HMM: window=" + winSize);
		HMM hmm = new HMM();
		hmm.train(trainAlphas, trainLabel);
		hmm.test(testAlphas, testLabel);
	}
	
	/**
	 * Reads aligned sequences from file.
	 * @return map [animal]->[sequence] 
	 */
	public static TreeMap<String, String> readSequences(String filePath){
		TreeMap<String, String> result = new TreeMap<>();
		try {
			Scanner in = new Scanner(new File(filePath));
			String animal = in.next().substring(1);
			StringBuilder seq = new StringBuilder();
			while(in.hasNext()){
				String line = in.next();
				if(line.startsWith(">")){
					result.put(animal, seq.toString());
					animal = line.substring(1);
					seq = new StringBuilder();
				}
				else{
					seq.append(line);
				}
			}
			result.put(animal, seq.toString());
			in.close();
		}
		catch (FileNotFoundException e) { e.printStackTrace(); }

		return result;
	}

}