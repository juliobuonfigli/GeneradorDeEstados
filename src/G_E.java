
import ij.plugin.PlugIn;
import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.measure.ResultsTable;
import java.text.DecimalFormat;
import java.awt.*;
import ij.plugin.*;
import ij.plugin.frame.*;
import java.util.Random;
import java.math.BigDecimal;
import java.math.BigInteger;

public class G_E implements PlugIn {
    
	private String[] ssd1={"NS", "SE", "SC"};
    private String[] ssd2={"NS", "SE", "SC", "Any"};
    private String[] ssd3={"SC-SC-SC", "NS-NS-NS", "SE-SE-SE", "SC-NS-SE", "SC-SC-NS", "SC-SC-SE",
    		              "SC-SE-SE", "SC-NS-NS", "SE-SE-NS", "SE-NS-NS", "Any"};
    private String[] sl={"0.1", "0.05", "0.01"};
    private int w, h, max, iter, seed;
    private String Srg, Sd2, Sd3, St, Sg, Ssl;
    private double SL;
    private boolean[] r, g, b;
    private Random rand;
    private ImagePlus img;
    private static DecimalFormat df1, df2;
    //private BigInteger fact; 
    

public void mainDialog()
	{
	GenericDialog gd = new GenericDialog("Multiple Relationship Finder");
	gd.addChoice("RedGreen: ", ssd1, ssd1[0]);
	gd.addChoice("Double 2: ",  ssd1, ssd1[0]);
	gd.addChoice("Double 3: ",  ssd1, ssd1[0]);
	gd.addChoice("Triples: ", ssd3, ssd3[0]);
	gd.addChoice("General: ",  ssd2, ssd2[0]);
	gd.addChoice("Significance level: ", sl, sl[0]);
	gd.addNumericField("Width: ", 10, 1);
	gd.addNumericField("Height: ", 10, 1);
	gd.addNumericField("Maximum: ", 60, 1);
	gd.addNumericField("Iteration (x1000): ", 100, 4);
	gd.addNumericField("Seed: ", 1, 9);
	gd.showDialog();
	if(gd.wasCanceled()) 
		return;
	Srg = gd.getNextChoice();
    Sd2 = gd.getNextChoice();
    Sd3 = gd.getNextChoice();
    St = gd.getNextChoice();
    Sg = gd.getNextChoice();
    Ssl = gd.getNextChoice();
    w = (int)gd.getNextNumber();
    h = (int)gd.getNextNumber();
    max = (int)gd.getNextNumber();
    iter = (int)gd.getNextNumber()*1000;
    seed = (int)gd.getNextNumber();
    }

public ImagePlus NEWFIGURE(boolean[] vector, String nombre)
	{
	ImagePlus img = IJ.createImage(nombre, "8-bit Black", w, h, 1); 	
	ImageProcessor IMG = img.getProcessor();
	int i=0;
	for(int y=0; y<w; y++)
		{
		for(int x=0; x<h; x++)
			{
			if(vector[i]==true)
				IMG.putPixelValue(x, y, 255); 
			i++;
			}
		} 
	img.setProcessor(IMG);
	return img;
	}

public BigInteger FAC(int s)
	{
	BigInteger x=BigInteger.valueOf(s);
	BigInteger y;
	if(s<=1) 
		x=BigInteger.valueOf(1); 
	else 
		{
		for(int i=s-1; i>0; i--)
			{
			y=BigInteger.valueOf(i);
			x=x.multiply(y); 
			}
		}
	return x;
	}

public double BinCoef(int n, int k)
	{
	BigInteger den=FAC(k).multiply(FAC(n-k));
	double res=(FAC(n).divide(den)).doubleValue();
	return res;	
	}

public double HG(int d, int n, int N, int x)  //funcion hipergeometrica
	{
	double res=(double)(BinCoef(d, x)*BinCoef(N-d, n-x)/BinCoef(N, n));
	return res;
	}

public double aHG(int d, int n, int N, int x)  //funcion hipergeometrica acumulada
	{
	double res=0;
	for(int i=x; i<n+1; i++)
		res=res+HG(d, n, N, i);
	return res; 	
	}

public double EM(int d, int n, int N)  //Esperanza matematica
{
double res=0;
for(int i=1; i<n+1; i++)
	res=res+i*HG(d, n, N, i);
return res; 	
}

public int DOBLE(boolean[] c1, boolean[] c2)
	{
	int n=0;
	for(int i=0; i<w*h; i++)
		{
		if(c1[i]==true && c2[i]==true)
			n++;
		}
	return n;
	}

public int TRIPLE(boolean[] c1, boolean[] c2, boolean[] c3)
	{
	int n=0;
	for(int i=0; i<w*h; i++)
		{
		if(c1[i]==true && c2[i]==true && c3[i]==true)
			n++;
		}
	return n;
	}

public boolean[] LOAD(boolean[] vec, int n, Random r)
	{
	int pos;
	for(int i=0; i<w*h; i++)
		vec[i]=false;
	for(int j=0; j<n; j++)
		{
		do  {
			pos=Math.abs(r.nextInt(w*h));	
			}
			while(vec[pos]==true);
		vec[pos]=true;
		}
	return vec;	
	}

public int NumLet(String let)
	{
	int res;
	if(let.equals("SC"))
		res=2;
	else
		{
		if(let.equals("NS"))
			res=0;
		else
			res=-3;
		}	
	return res;
	}

public int NumNum(double s, double e, double sl)
	{
	int res;
	if(s<sl)
		res=2;
	else
		{
		if(e<sl)
			res=-3;
		else
			res=0;
		}		
	return res;
	}

public String LetNum(int num)
	{
	String res;
	if(num==2)
		res="SC";
	else
		{
		if(num==0)
			res="NS";
		else
			res="SE";
		}	
	return res;
	}

public void createRT(double[] c, double[] sc, double[] se, String[] ss, int it, int sem)
	{
	df2 = new DecimalFormat("0.00E0");
	df1 = new DecimalFormat("0.00");
	ResultsTable rt=new ResultsTable();
	rt.incrementCounter();
	rt.addValue("Results", "Coefficients");
	rt.addValue("RedGreen", df1.format(c[0]));
	rt.addValue("GreenRed", df1.format(c[1]));
	rt.addValue("RedBlue", df1.format(c[2]));
	rt.addValue("BlueRed", df1.format(c[3]));
	rt.addValue("GreenBlue", df1.format(c[4]));
	rt.addValue("BlueGreen", df1.format(c[5]));
	rt.addValue("TRed", df1.format(c[6]));
	rt.addValue("TGreen", df1.format(c[7]));
	rt.addValue("TBlue", df1.format(c[8]));
	rt.addValue("TGeneral", df1.format(c[9]));
	rt.addValue("Iterations", it);
	rt.incrementCounter();
	rt.addValue("Results", "Significances");
	rt.addValue("RedGreen", ss[0]);
	rt.addValue("GreenRed", "  ");
	rt.addValue("RedBlue", ss[1]);
	rt.addValue("BlueRed", " ");
	rt.addValue("GreenBlue", ss[2]);
	rt.addValue("BlueGreen", " ");
	rt.addValue("TRed", ss[3]);
	rt.addValue("TGreen", ss[4]);
	rt.addValue("TBlue", ss[5]);
	rt.addValue("TGeneral", ss[6]);
	rt.addValue("Iterations", sem);
	rt.incrementCounter();
	rt.addValue("Results", "pSC");
	rt.addValue("RedGreen", df2.format(sc[0]));
	rt.addValue("GreenRed", " ");
	rt.addValue("RedBlue", df2.format(sc[1]));
	rt.addValue("BlueRed", " ");
	rt.addValue("GreenBlue", df2.format(sc[2]));
	rt.addValue("BlueGreen", " ");
	rt.addValue("TRed", df2.format(sc[3]));
	rt.addValue("TGreen", df2.format(sc[4]));
	rt.addValue("TBlue", df2.format(sc[5]));
	rt.addValue("TGeneral", df2.format(sc[6]));
	rt.addValue("Iterations", " ");
	rt.incrementCounter();
	rt.addValue("Results", "pSE");
	rt.addValue("RedGreen", df2.format(se[0]));
	rt.addValue("GreenRed", " ");
	rt.addValue("RedBlue", df2.format(se[1]));
	rt.addValue("BlueRed", " ");
	rt.addValue("GreenBlue", df2.format(se[2]));
	rt.addValue("BlueGreen", " ");
	rt.addValue("TRed", df2.format(se[3]));
	rt.addValue("TGreen", df2.format(se[4]));
	rt.addValue("TBlue", df2.format(se[5]));
	rt.addValue("TGeneral", df2.format(se[6]));
	rt.addValue("Iterations", " ");
	rt.showRowNumbers(false);
	rt.show("Results"); 
	}
	
@Override
public void run(String arg0)
	{	
	mainDialog();
	
	SL = Double.parseDouble(Ssl);
	r=new boolean[w*h];
	g=new boolean[w*h];
	b=new boolean[w*h];

	rand=new Random(seed); 
		
	int[] DOBLES=new int[3];
	
	DOBLES[0]=NumLet(Srg);
	DOBLES[0]=NumLet(Sd2);
	DOBLES[0]=NumLet(Sd3);
	
	int cond1=DOBLES[0]+DOBLES[1]+DOBLES[2];
	int cond2=0;
	int cond3=0;
	
	if(St.equals("SC-SC-SC")) cond2=6;
	if(St.equals("NS-NS-NS")) cond2=0;
	if(St.equals("SE-SE-SE")) cond2=-9;
	if(St.equals("SC-NS-SE")) cond2=-1;
	if(St.equals("SC-SC-NS")) cond2=4;
	if(St.equals("SC-SC-SE")) cond2=1;
	if(St.equals("SC-SE-SE")) cond2=-4;
	if(St.equals("SC-NS-NS")) cond2=2;
	if(St.equals("SE-SE-NS")) cond2=6;
	if(St.equals("SE-NS-NS")) cond2=-3;
	if(St.equals("Any")) cond2=10;

	if(Sg.equals("SC")) cond3=2;
	if(Sg.equals("NS")) cond3=0;
	if(Sg.equals("SE")) cond3=-3;
	if(Sg.equals("Any")) cond3=10;

	double[] SS=new double[7];
	double[] eSS=new double[7];
	int[] ss=new int[7];
	int total=0;
	int nr, ng, nb, rg, rb, gb, rgb, dd, COND0, COND1, COND2, co1, co2, co3, condition, cont;
	nr=0; ng=0; rg=0;	
	
	nr=Math.abs(rand.nextInt(max))+1;
	ng=Math.abs(rand.nextInt(max))+1;
	r=LOAD(r, nr, rand);
	g=LOAD(g, ng, rand);
	rg=DOBLE(r, g);
	SS[0]=aHG(ng, nr, w*h, rg);
	eSS[0]=aHG(w*h-ng, nr, w*h, nr-rg);
	COND0=NumNum(SS[0], eSS[0], SL);
	
	/*System.out.println("nr: "+nr);
	System.out.println("ng: "+ng);
	System.out.println("w*h: "+w*h);
	System.out.println("rg: "+rg);
	System.out.println("SS[0]: "+SS[0]);
	System.out.println("eSS[0]: "+eSS[0]);
	System.out.println("HG: "+HG(ng, nr, w*h, rg));
	System.out.println("FAC: "+FAC(29));
	System.out.println("BinCoef: "+BinCoef(nr, rg));*/
	do 
		{	
		cont=0;
		condition=0;  
		if(total%100==0)
			{
			do
				{
				nr=rand.nextInt(max)+1;
				ng=rand.nextInt(max)+1;
				r=LOAD(r, nr, rand);
				g=LOAD(g, ng, rand);
				rg=DOBLE(r, g);
				SS[0]=aHG(ng, nr, w*h, rg);
				eSS[0]=aHG(w*h-ng, nr, w*h, nr-rg);
				COND0=NumNum(SS[0], eSS[0], SL);
				cont++;
				} 
			while(COND0!=DOBLES[0] && cont<100);
			}
		nb=rand.nextInt(max)+1;
		b=LOAD(b, nb, rand);
		//	aHG(N-d, N-n, N, n-x)
		rb=DOBLE(r, b);
		gb=DOBLE(g, b);
		rgb=TRIPLE(r, g, b);
		//aHG(d, n, N, x) 
		SS[1]=aHG(nb, nr, w*h, rb);
		SS[2]=aHG(ng, nb, w*h, gb);
		eSS[1]=aHG(w*h-nb, nr, w*h, nr-rb);
		eSS[2]=aHG(w*h-ng, nb, w*h, nb-gb);
		
		SS[3]=aHG(gb, nr, w*h, rgb);
		SS[4]=aHG(rb, ng, w*h, rgb);
		SS[5]=aHG(rg, nb, w*h, rgb);
		eSS[3]=aHG(w*h-gb, nr, w*h, nr-rgb);
		eSS[4]=aHG(w*h-rb, ng, w*h, ng-rgb);
		eSS[5]=aHG(w*h-rg, nb, w*h, nb-rgb);
	
		dd=(int)Math.round(nr*ng/(w*h));
				
		SS[6]=aHG(dd, nb, w*h, rgb);
		eSS[6]=aHG(w*h-dd, nb, w*h, nb-rgb);
	
		for(int i=0; i<7; i++)
			ss[i]=NumNum(SS[i], eSS[i], SL);
			
		COND1=ss[0]+ss[1]+ss[2];
		COND2=ss[3]+ss[4]+ss[5];
		
		co1=0; co2=0; co3=0;
		if(cond1==COND1)
			co1=1;
		if(cond2==COND2 || cond2==10)
			co2=1;
		if(cond3==ss[6] || cond3==10)
			co3=1;
		
		condition=co1+co2+co3;
		
	   	total++;
		}
		while(condition<3 && total<iter);
	
	ImagePlus RED=NEWFIGURE(r, "Red");	RED.show();
	ImagePlus GREEN=NEWFIGURE(g, "Green"); GREEN.show();
	ImagePlus BLUE=NEWFIGURE(b, "Blue"); BLUE.show();
	
	IJ.run(img, "Merge Channels...", "c1=Red c2=Green c3=Blue create");
	IJ.run(img, "Size...", "width=80 height=80 constrain average interpolation=None");
	//img.show();
	
	String[] sS=new String[7];
	for(int i=0; i<7; i++)
		sS[i]=LetNum(ss[i]);
	    	
	double[] c=new double[10];
	c[0]=(double)rg/nr; c[1]=(double)rg/ng; c[2]=(double)rb/nr; c[3]=(double)rb/nb; c[4]=(double)gb/ng; 
	c[5]=(double)gb/nb; c[6]=(double)rgb/nr; c[7]=(double)rgb/ng; c[8]=(double)rgb/nb; c[9]=(double)3*rgb/(nr+ng+nb); 
	createRT(c, SS, eSS, sS, total, seed);
	
	/*for(int i=0; i<7; i++) {
		System.out.println("e: "+eSS[i]);
		System.out.println("c: "+SS[i]); }*/		
	}

}

/*
public static void main(String[] args) 
	{
	// set the plugins.dir property to make the plugin appear in the Plugins menu
	Class<?> clazz = G_E.class;
	String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
	String pluginsDir = url.substring("file:".length(), url.length() - clazz.getName().length() - ".class".length());
	System.setProperty("plugins.dir", pluginsDir);
	// start ImageJ
	new ImageJ();
	// run the plugin
	IJ.runPlugIn(clazz.getName(), "");
	}

}


/* factorial recursivo
 * public long FAC(int s)
	{
    if(s<=1) 
        return 1;
    else 
        return s*FAC(s-1);
    }
 */

