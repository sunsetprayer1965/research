package computewritec;

public class test {
    public static void main(String[] args) throws Exception {
        // First compare the results using 1.txt and 2.txt from Panel_1.
        double[][] v1 = ComputeWriteC3.loadChromosomeFromFile(
            "C:\\Users\\weilu\\Downloads\\Panels (1)\\Panels\\Panel_1\\1.txt");
        double[][] v2 = ComputeWriteC3.loadChromosomeFromFile(
            "C:\\Users\\weilu\\Downloads\\Panels (1)\\Panels\\Panel_1\\2.txt");

        int K = 300; // Number of random samples to compare.
        double[][][] pts = ComputeWriteC3.pregenRandomPoints(20241011L, K);

        int mismatch = 0;
        for (int k = 0; k < K; k++) {
            double[] p1 = pts[k][0], p2 = pts[k][1];

            double[][] v1Old = ComputeWriteC3.addedVector(v1, p1);
            double[][] v2Old = ComputeWriteC3.addedVector(v2, p2);
            double old01 = ComputeWriteC3.compute_link_state(v1Old, v2Old);

            double new01 = ComputeWriteC3.compute_link_state_augmented(v1, v2, p1, p2);

            System.out.printf("k=%d old=%.1f new=%.1f%n", k, old01, new01);

            if (old01 != new01) {
                System.out.println("⚠️  Mismatch at k=" + k);
                mismatch++;
            }
        }
        System.out.println("✅ 对拍完成，mismatch=" + mismatch);
    }
}
