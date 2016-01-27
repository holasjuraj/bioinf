public class LogNum {
	
	/** Value of log(x) */
	private double logx = 0;
	private boolean isZero = false;
	
	public LogNum(double x){
		if(x <= 0){
			isZero = true;
		}else{
			logx = Math.log(x);
		}
	}
	public LogNum(){
		this(0);
	}
	
	/** Multiplication */
	public static LogNum mul(LogNum a, LogNum b){
		LogNum result = new LogNum();
		result.isZero = a.isZero || b.isZero;
		result.logx = a.logx + b.logx;
		return result;
	}
	public static LogNum mul(double a, LogNum b){
		return mul(new LogNum(a), b);
	}
	public static LogNum mul(LogNum a, LogNum b, LogNum c){
		return mul(mul(a, b), c);
	}
	
	/** Comparison of two numbers */
	public boolean isGt(LogNum b){
		return !isZero && (b.isZero || logx > b.logx);
	}
	
	public String toString(){
		if(isZero){ return "0.0"; }
		// From JavaDoc: for values of x near 0, the exact sum of expm1(x)+1 is
		// much closer to the true result of e^x than exp(x).
		return new Double(Math.expm1(logx)+1).toString();
	}

}