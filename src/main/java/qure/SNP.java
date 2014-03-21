package qure;

public class SNP
{
	public double position;
	public char consensus;
	public char base;
	
	public SNP()
	{
		position=-1;
		consensus='X';
		base='X';
	}
	
	public SNP(double p, char c, char b)
	{
		position=p;
		consensus=c;
		base=b;
	}
	
	public String toString()
	{
		String s = consensus+"_"+position+"_"+base+",";
		return s;
	}
	
	public boolean equals(SNP b)
	{
		return (this.position==b.position && this.consensus==b.consensus && this.base==b.base);
	}
}
