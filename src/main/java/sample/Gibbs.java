package sample;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.Random;

/**
 *
 * Gibb's Sampling Steps:<br>
 * 1. Set every variable to a random value.<br>
 * 2. Choose a variable to update. <br>
 * 3. Randomly Select (aka "Sample") a new value for the variable based on the
 * current conditions. <br>
 * 4. Repeat from Step 2.
 *
 *
 */
public class Gibbs {

    private static ArrayList<String> sequences = new ArrayList<String>();
    Hashtable<String, Integer> start;
    private int motifLength;

    /**
     * Constructs and performs Gibb's Sampling in order to find repeated motifs.
     *
     * @param seq
     *            A String array of the sequences that will be used.
     * @param motifLength
     *            An Integer that shows the length of the motif or pattern we
     *            are trying to find, this value is given.
     */
    public Gibbs(String[] seqArray, int motifLength) {
        sequences.addAll(Arrays.asList(seqArray));
        this.motifLength = motifLength;
        this.start = generateRandomValue();
        sample();
        System.out.println(start);
    }

    /**
     * This method is repeated 2000 times.
     *
     * @param start
     *            A HashTable containing the sequence as a key, and the random
     *            integer to be used as the value.
     */
    private void sample() {
        for (int j = 0; j < 2000; j++) {
            Random rand = new Random();
            int chosenSeqIndex = rand.nextInt(sequences.size());
            String chosenSequence = sequences.get(chosenSeqIndex);
            ArrayList<Double> scores = new ArrayList<Double>();
            // i = possibleStart
            for (int i = 0; i < chosenSequence.length() - motifLength + 1; i++) {
                String tempMotif = chosenSequence.substring(i, i + motifLength);
                double p = calculateP(tempMotif, chosenSeqIndex);
                double q = calculateQ(tempMotif, chosenSeqIndex, i);
                scores.add(q / p);
            }
            double sum = 0;
            for (double d : scores) {
                sum += d;
            }
            for (int i = 0; i < scores.size(); i++) {
                scores.set(0, scores.get(i) / sum);
            }

            double random = rand.nextDouble();
            double dubsum = 0;
            for (double d : scores) {
                dubsum += d;
                if (random == dubsum) {
                    start.put(chosenSequence, scores.indexOf(d));
                }
            }
        }

    }

    /**
     * Calculates the probability of a letter in this position.
     *
     * @param tempMotif
     *            The motif being used for this calculation.
     * @param chosenSeqIndex
     *            The index of the sequence being used for this calculation,
     *            useful for skipping all of this sequences calculations and
     *            focusing on the other ones.
     * @return A double of the probability of a letter in this position.
     */
    private double calculateQ(String tempMotif, int chosenSeqIndex,
                              int possibleStart) {
        double q = 1;
        int start = possibleStart;
        int end = possibleStart + tempMotif.length();
        double denominator = sequences.size() - 1;
        for (String s : sequences) {
            double numerator = 0;
            if (s.equals(sequences.get(chosenSeqIndex)))
                continue;
            if (end > s.length()) {
                q *= 0.01;
                continue;
            }
            String thisMotif = s.substring(start, end);
            char[] letters = tempMotif.toCharArray();
            char[] seqLetters = thisMotif.toCharArray();

            for (int i = 0; i < tempMotif.length(); i++) {
                if (letters[i] == seqLetters[i])
                    numerator++;
            }
            if (numerator == 0)
                q *= 0.01;
            else
                q *= (numerator / denominator);
        }
        return q;
    }

    /**
     * Calculates the probability of a letter randomly selected.
     *
     * To find this value, the method loops through each letter of the selected
     * temporary motif, and loops through the other sequences. While looping
     * through the other sequences, we find the amount of same letters in each
     * other sequence, along with the total length of all other sequences. The
     * value P is a product of every result, each result being the amount of
     * letters of the same kind over the total amount of letters.
     *
     * @param tempMotif
     *            The motif being used for this calculation.
     * @param chosenSeqIndex
     *            The index of the sequence being used for this calculation,
     *            useful for skipping all of this sequences calculations and
     *            focusing on the other ones.
     * @return A double of the probability of a letter randomly selected.
     */
    private double calculateP(String tempMotif, int chosenSeqIndex) {
        double p = 1;
        for (char c : tempMotif.toCharArray()) {
            double sameLetters = 0;
            double totalLength = 0;
            for (String s : sequences) {
                if (s.equals(sequences.get(chosenSeqIndex)))
                    continue;
                char[] seqLetters = s.toCharArray();
                for (char x : seqLetters)
                    if (c == x)
                        sameLetters++;
                totalLength += s.length();
            }
            p *= (sameLetters / totalLength);
        }
        return p;
    }

    /**
     * Calculates and stores every random value. Generates a random from 0 to a
     * value of each individual sequences length subtracted by the motif length.
     *
     * @return A HashTable containing the sequence as a key, and the random
     *         integer to be used as the value.
     */
    private Hashtable<String, Integer> generateRandomValue() {
        Random rand = new Random();
        Hashtable<String, Integer> randomValues = new Hashtable<String, Integer>();
        for (String seq : sequences) {
            int randomVal = rand.nextInt(seq.length() - motifLength);
            randomValues.put(seq, randomVal);
        }
        return randomValues;
    }

    public static void main(String[] args) {
        String[] data = { "ABCDAAAABDB", "AAAADCBBCA", "DDBCABAAAACBBD",
                "AABAAAACCDD" };
        int length = 4;
        @SuppressWarnings("unused")
        Gibbs gibbs = new Gibbs(data, length);
    }

}
