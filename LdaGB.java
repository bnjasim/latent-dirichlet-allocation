/*
 * LdaGB is free software; you can redistribute it and/or modify it
 * Created on Jul 16, 2013
 * @author Binu Jasim: bnjasim@gmail.com
 * Vocabulary size should be <= 32,767
 */
import java.util.Scanner;
import java.io.*;

/**
 * The Gibbs sampling implementation of Latent Dirichlet Allocation
 *
 * @author binu
 */
public class LdaGB {

    private final int MAX_ITER = 100;
    private final double alpha = 0.5; //Uniform prior
    private final double beta = 0.1;

    private final int D; // Number of Documents
    private final int V; // Size of Vocabulary
    private final int Z; // Number of Topics
    private final int W; // Total number of words/terms in the whole corpus

    private final short ArrayD[]; // Document each word belongs to - size = W
    private final short ArrayV[]; // Vocabulary each word corresponds to - size = W
    private short ArrayZ[]; // Topic each word corresponds to - size = W

    private int TableDZ[][]; // size DxZ
    private int TableZV[][]; // size ZxV
    private int RowSumZV[];

    /* Initialization of the Tables */
    LdaGB(String file, int d, int v, int z)
    {
        Scanner in = null;
        try {
            in = new Scanner(new File(file));
        }
        catch(Exception e){
            System.out.println("The Input document Doesn't exist");
            System.exit(1);
        }
        D = d;
        V = v;
        Z = z;

        // First estimate V by readin the doc file
        int wc = 0;
        Scanner lScan = null;
        for(int i=0;i<D;i++)
        {
            String line = in.nextLine();
            lScan = new Scanner(line);
            while(lScan.hasNext())
            {
                String s = lScan.next();
                String[] t = s.split(":");
                //System.out.println(t[0]);
                wc += Integer.parseInt(t[1]);
            }
        }
        lScan.close();
        in.close();
        W = wc;
        System.out.println("W = "+W);
        ArrayD = new short[W];
        ArrayV = new short[W];
        ArrayZ = new short[W];
        TableDZ = new int[D][Z];
        TableZV = new int[Z][V];
        RowSumZV = new int[Z];
        // Now Initialize the Tables and Arrays
        // Reading the doc file for one last time
        wc = 0;
        try {
            in = new Scanner(new File(file));
        }
        catch(Exception e){
            System.out.println("Something wrong with Reading document");
            System.exit(1);
        }
        for(int i=0;i<D;i++)
        {
            String line = in.nextLine();
            lScan = new Scanner(line);
            while(lScan.hasNext())
            {
                String s = lScan.next();
                String[] t = s.split(":");
                int f = Integer.parseInt(t[1]);
                for(int j=0; j<f; j++)
                {
                    z = sampleTopicUniform();
                    short voc = (short) (Integer.parseInt(t[0]) - 1);
                    ArrayV[wc] = voc;
                    ArrayZ[wc] = (short)z;
                    ArrayD[wc] = (short)i;
                    TableDZ[i][z]++;
                    TableZV[z][voc]++;
                    RowSumZV[z]++;
                    wc++;
                }
            }
        }
        lScan.close();
        in.close(); 
        /*
        for(int i=0;i<wc;i++)
        {
            System.out.println(ArrayD[i]+" "+ArrayV[i]+" "+ArrayZ[i]);
        }
        */
        
    } // End of the Constructor



    private void gibbsIterate() {

        double p[] = new double[Z];

        for(int i=0; i<MAX_ITER; i++)
        {
            for(int wc=0; wc<W; wc++)
            {
                int curr_topic = ArrayZ[wc];
                int curr_voc = ArrayV[wc];
                int curr_doc = ArrayD[wc];

                for(int k=0; k<Z; k++)
                    p[k] = (TableDZ[curr_doc][k] - (curr_topic==k?1:0) + alpha) * (TableZV[k][curr_voc] - (curr_topic==k?1:0) + beta) / (RowSumZV[k] - (curr_topic==k?1:0) + beta*V);
                    
                short z = sampleTopicP(p);
                if(z != curr_topic) // have to alter the tables
                    
                {
                    ArrayZ[wc] = z;
                    TableDZ[curr_doc][curr_topic]--;
                    TableDZ[curr_doc][z]++;
                    TableZV[curr_topic][curr_voc]--;
                    TableZV[z][curr_voc]++;
                    RowSumZV[curr_topic]--;
                    RowSumZV[z]++;
                }
            }
            System.out.println("Iteration = "+(i+1));
        }
    }

    private short sampleTopicP(double[] p) {
        short z = 0;
        double r = Math.random();

        // Normalize p
        double sum = 0;
        double cum[] = new double[Z];

        for(int i=0; i<Z; i++)
        {
            sum += p[i];
        }
        p[0] = p[0]/sum;
        cum[0] = p[0];
        for(int i=1; i<Z; i++)
        {
            p[i] = p[i]/sum;
            cum[i] = p[i] + cum[i-1];
        }
/* for(int zz= 0; zz<Z; zz++){
        System.out.print(cum[zz] + " "); } System.out.println(); */
        for(short i=0; i<Z; i++)
        {
            if(r < cum[i])
            {
                z = i;
                break;
            }
        }
        //System.out.println(z);
        return z;
    }

    private int sampleTopicUniform() {
        int z = 0;
        double r = Math.random();

        for(int i=0; i<Z; i++ )
        {
            if(r < (i+1)/(double)Z)
            {
                z = i;
                break;
            }
        }
        //System.out.println(z);
        return z;
    }

    private void printPhi()
    {
        try {
            PrintWriter out = new PrintWriter(new File("../outPhi.txt"));

            for(int i=0; i<Z; i++)
            {
                for(int j=0; j<V; j++)
                {
                    double b = (TableZV[i][j] + beta) / (RowSumZV[i] + beta*V);
                    String s = String.format("%f ", b);
                    out.write(s);
                }
                out.write("\n");
            }
            out.close();
        }
        catch(Exception e) {
            System.out.println("Couldn't save the result");
        }
    }

    private void check() {
        for(int i=0;i<Z; i++)
        {
            System.out.println(RowSumZV[i]);
        }

        //check whether any negative in TableZV
        int c = 0;
        for(int i=0; i<Z; i++)
        {
            for(int j=0; j<V; j++)
            {
                System.out.printf("%d ",TableZV[i][j]);
               
            }
            System.out.println();
        }
    }

    public static void main(String[] args)
    {
        String fileDoc = "../corpus.txt";
        int d = 100;
        int v = 25;
        int z = 8;

        LdaGB lda = new LdaGB(fileDoc, d, v, z);
        lda.gibbsIterate();
        //lda.check();
        lda.printPhi();
    }


}
