package sample;

import java.io.*;
import java.util.*;

public class MetropolisHastings
{
    public static void main(String[] args) throws IOException
    {
        final int NUM_ITERATIONS = 1000000;

        File file = new File("data/data.txt");
        Scanner scan = new Scanner(file);

        int dim1 = scan.nextInt();
        int dim2 = scan.nextInt();
        scan.nextLine();

        HashMap<Double, Integer> chiSquareValues = new HashMap<Double, Integer>();

        int[] observed = new int[dim1 * dim2];

        for(int i = 0; i < dim1 * dim2; i++)
        {
            observed[i] = scan.nextInt();
        }

        int[] columnSums = new int[dim1];
        int total = 0;

        //add up columns
        for(int i = 0; i < dim1; i++)
        {
            for(int j = 0; j < dim2; j++)
            {
                total += observed[i * dim2 + j];
            }

            columnSums[i] = total;
            total = 0;
        }

        int[] rowSums = new int[dim2];

        //add up rows
        for(int i = 0; i < dim2; i++)
        {
            for(int j = 0; j < dim1; j++)
            {
                total += observed[j * dim2 + i];
            }

            rowSums[i] = total;
            total = 0;
        }

        int n = 0;

        //get total
        for(int i = 0; i < dim2; i++)
        {
            n += columnSums[i];
        }

        //null hypothesis: expected
        double[] expected = new double[dim1 * dim2];

        //calculate expected independence values
        for(int i = 0; i < dim1; i++)
        {
            for(int j = 0; j < dim2; j++)
            {
                expected[i * dim1 + j] = (columnSums[i] * rowSums[j]) / (double)(n);
            }
        }

        //calc Chi-Square for observed compared to expected
        double obs = calcChiSquare(observed, expected);

        int dim3 = scan.nextInt();
        int dim4 = scan.nextInt();
        scan.nextLine();

        int[][] markovBasis = new int[dim3][dim4];

        //store moves of markov basis
        for(int i = 0; i < dim3; i++)
        {
            for(int j = 0; j < dim4; j++)
            {
                markovBasis[i][j] = scan.nextInt();
            }

            if(scan.hasNextLine())
                scan.nextLine();
        }

        //algorithm

        Random rand = new Random();
        int markovPick = 0;
        int[] tempMarkov = new int[dim4];
        int epsPick = 0;
        boolean inFiber = true;
        double u = 0;
        int sig = 0;

        //copy observed to x
        int[] x = new int[observed.length];

        for(int i = 0; i < observed.length; i++)
        {
            x[i] = observed[i];
        }

        for(int i = 0; i < NUM_ITERATIONS; i++)
        {
            //initialize
            inFiber = true;

            //step 2
            markovPick = rand.nextInt(dim3);
            epsPick = rand.nextInt(2);

            for(int j = 0; j < dim4; j++)
            {
                tempMarkov[j] = (int)(Math.pow(-1, epsPick) * markovBasis[markovPick][j]);
            }

            //step 3: test for any negative values
            for(int j = 0; j < x.length; j++)
            {
                if(x[j] + tempMarkov[j] < 0)
                    inFiber = false;
            }

            //step 4
            if(inFiber)
            {
                u = rand.nextDouble();

                if(u < calcTransProb(x, tempMarkov))
                {
                    for(int k = 0; k < x.length; k++)
                    {
                        x[k] = x[k] + tempMarkov[k];
                    }
                }
            }

            double chiSquare = calcChiSquare(x, expected);

            //step 5
            if(chiSquare >= obs)
                sig++;
        }

        System.out.println("p-value: " + (double)sig / NUM_ITERATIONS);
    }

    public static double calcTransProb(int[] x, int[] markov)
    {
        double prob = 1;

        for(int i = 0; i < x.length; i++)
        {
            if(markov[i] > 0)
            {
                for(int j = 1; j <= markov[i]; j++)
                {
                    prob = prob / (x[i] + j);
                }
            }

            if(markov[i] < 0)
            {
                for(int j = 0; j < markov[i] * -1; j++)
                {
                    prob = prob * (x[i] - j);
                }
            }
        }

        return prob;
    }

    public static double calcChiSquare(int[] observed, double[] expected)
    {
        double total = 0;

        for(int i = 0; i < observed.length; i++)
        {
            total += (Math.pow(((double)observed[i] - expected[i]), 2)) / expected[i];
        }

        return total;
    }
}
