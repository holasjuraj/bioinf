import java.util.List;

public class HMM {

	final static char	state1 = 'n',	// State "not exon"
						state2 = 'e';	// State "exon"
	LogNum[] pi;	// Initial probabilities
	LogNum[][] t;	// Transition probabilities
	LogNum[][] e;	// Emission probabilities

	/**
	 * Training of model: counting all probabilities
	 * @param alphas : training set of alphas (emissions)
	 * @param states : training set of states
	 */
	public void train(List<Integer> alphas, String states){
		t = new LogNum[2][2];
		e = new LogNum[2][21];
		pi = new LogNum[2];
		
		double  s1count = 0, s2count = 0,   // States count
				change12 = 0, change21 = 0; // Changes count
		double[] a1count = new double[21],  // Count alphas in state 1
				 a2count = new double[21];  // Count alphas in state 2
		char lastState = state1;
		
		// Count emissions, states and state changes
		for(int i = 0; i < alphas.size(); i++){
			int a = alphas.get(i);
			char s = states.charAt(i);
			
			if(s == state1){
				s1count++;
				if(lastState == state2){ change21++; }
				a1count[a]++;
			}
			else{
				s2count++;
				if(lastState == state1){ change12++; }
				a2count[a]++;
			}
			lastState = s;
		}
		
		// Initial probabilities
		pi[0] = new LogNum(s1count / (s1count + s2count));
		pi[1] = new LogNum(s2count / (s1count + s2count));
		
		// Transition probabilities
		t[0][0] = new LogNum((s1count - change12) / s1count);
		t[0][1] = new LogNum(change12 / s1count);
		t[1][1] = new LogNum((s2count - change21) / s2count);
		t[1][0] = new LogNum(change21 / s2count);
		
		// Emission probabilities
		for(int i = 0; i < 21; i++){
			e[0][i] = new LogNum(a1count[i] / s1count);
			e[1][i] = new LogNum(a2count[i] / s2count);
		}

		// Print
		System.out.println("Initial probabilities:");
		System.out.println(pi[0]);
		System.out.println(pi[1]);
		System.out.println("Transition probabilities:");
		System.out.println(t[0][0] +"\t"+ t[0][1]);
		System.out.println(t[1][0] +"\t"+ t[1][1]);
		System.out.println("Emission probabilities:");
		for(int i = 0; i < 21; i++){
			System.out.println(e[0][i] +"\t"+ e[1][i]);
		}
		System.out.println();
	}
	
	/**
	 * Implementation of Viterbi algorithm, using previously trained
	 * probabilities. Table of dynamic programming is composed of Fields,
	 * tracking both probabilities and paths through table.
	 * @param alphas : list of alphas (emissions)
	 * @return sequence of inferred states
	 */
	public String viterbi(List<Integer> alphas){
		/** Structure for dynamic programming */
		class field{
			LogNum p = new LogNum(); // Probability
			int from = 0;            // Previous field
			public field(){}
			public field(LogNum p, int from){
				this.p = p;
				this.from = from;
			}
		};
		
		// Viterbi algorithm
		field[][] A = new field[alphas.size()][2];
		A[0][0] = new field(LogNum.mul(pi[0], e[0][alphas.get(0)]), -1);
		A[0][1] = new field(LogNum.mul(pi[1], e[1][alphas.get(0)]), -1);
		
		for(int j = 1; j < alphas.size(); j++){
			int a = alphas.get(j);
			for(int k = 0; k < 2; k++){
				LogNum	p0 = LogNum.mul(A[j-1][0].p , t[0][k] , e[k][a]),
						p1 = LogNum.mul(A[j-1][1].p , t[1][k] , e[k][a]);
				A[j][k] = new field();
				if(p0.isGt(p1)){
					A[j][k].p = p0;
					A[j][k].from = 0;
				}
				else{
					A[j][k].p = p1;
					A[j][k].from = 1;
				}
			}
		}
		
		// Path reconstruction
		StringBuilder sb = new StringBuilder();
		int from = 0;
		if(A[A.length-1][1].p.isGt(A[A.length-1][0].p)){
			from = 1;
		}
		for(int i = A.length-1; i >= 0; i--){
			if(from == 0){
				sb.append(state1);
			}else{
				sb.append(state2);
			}
			from = A[i][from].from;
		}
		return sb.reverse().toString();		
	}
	
	/**
	 * Testing model on new data, computing its success rate.
	 * @param alphas : testing set of alphas (emissions) - known to model
	 * @param states : testing set of states - only to check results of model
	 */
	public void test(List<Integer> alphas, String states){
		// Launch Viterbi algorithm
		String guess = viterbi(alphas);		
		// Results
		int s1correct = 0, s1wrong = 0,
			s2correct = 0, s2wrong = 0;
		for(int i = 0; i < states.length(); i++){
			if(states.charAt(i) == state1){
				if(guess.charAt(i) == state1){ s1correct++; }
				else{ s1wrong++; }
			}
			if(states.charAt(i) == state2){
				if(guess.charAt(i) == state2){ s2correct++; }
				else{ s2wrong++; }
			}
		}
		System.out.println("Results:");
		System.out.println("state 1 correct = " + s1correct + "\tstate 1 wrong = " + s1wrong);
		System.out.println("state 2 correct = " + s2correct + "\tstate 2 wrong = " + s2wrong);
		System.out.println("Total correct = " + (double)(s1correct+s2correct)/(double)(states.length()) );
	}

}