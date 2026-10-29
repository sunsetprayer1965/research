package computewritec;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;

public class ComputeWriteC3 {

    // ====== 与 C2 一样的全局变量/常量 ======
    static int numberOfSample;
    static String[] listChromosomeFilePaths;

    // ✅ 你自己改成你想要的输出目录（保持 C2 风格：写死 Windows 路径）
    // static final String RESULT_DIR = "C:\\Users\\weilu\\Desktop\\research william luo\\python_script\\java\\RESULT";

    // ====== 路径构造（与 C2 一样：baseDir + "/1.txt" ...） ======
    public static String[] create_Strings(int simpleSize, String baseDir) {
        String[] solution = new String[simpleSize];
        for (int i = 0; i < simpleSize; i++) {
            int j = i + 1;
            solution[i] = baseDir + "/" + j + ".txt"; // 用 / 避免转义
        }
        return solution;
    }

    // ====== 统计该目录下连续存在的 1.txt, 2.txt, ...（与 C2 一样） ======
    public static int countSequentialTxtFiles(String baseDir) {
        int count = 0;
        while (true) {
            java.io.File f = new java.io.File(baseDir + "/" + (count + 1) + ".txt");
            if (f.exists() && f.isFile()) count++;
            else break;
        }
        return count;
    }

    // ====== main：结构/输出行为保持 C2 一样，仅把 Monte Carlo(k) 并行化 ======
    public static void main(String[] args) throws FileNotFoundException, IOException {


        if (args.length < 2) {
            System.out.println("Usage: java computewritec.ComputeWriteC3 <RESULT_DIR> <baseDir1> <baseDir2> ...");
            return;
        }

        String RESULT_DIR = args[0];

        String[] baseDirs = new String[args.length - 1];
        System.arraycopy(args, 1, baseDirs, 0, args.length - 1);

        int numberOfPoints = 1000;
        final int threads = 8;
        for (String baseDir : baseDirs) {

            // 1) 统计样本数
            numberOfSample = countSequentialTxtFiles(baseDir);
            if (numberOfSample == 0) {
                System.out.println("⚠️ 目录没有发现 1.txt, 2.txt... : " + baseDir);
                continue;
            }

            // 2) 构造文件列表
            listChromosomeFilePaths = create_Strings(numberOfSample, baseDir);

            // 3) label：目录最后一段
            String label = baseDir.replace('\\', '/');
            int lastSlash = label.lastIndexOf('/');
            label = (lastSlash >= 0) ? label.substring(lastSlash + 1) : label;

            // 4) 输出文件（每个目录一个 txt，和 C2 一样）
            java.io.File resultDir = new java.io.File(RESULT_DIR);
            resultDir.mkdirs();

            java.io.File outPath = new java.io.File(resultDir, label + "_linking_frequency.txt");

            System.out.println("\n========================");
            System.out.println("▶ Running on folder: " + baseDir);
            System.out.println("▶ Samples detected:  " + numberOfSample);
            System.out.println("▶ MonteCarlo N:      " + numberOfPoints);
            System.out.println("▶ Threads:           " + threads);
            System.out.println("▶ Output file:       " + outPath.getAbsolutePath());
            System.out.println("========================");

            // 5) 写文件：每行一个 pair 的百分比（与 C2 一样）
            try (FileWriter fw = new FileWriter(outPath);
                 PrintWriter out = new PrintWriter(new BufferedWriter(fw))) {

                for (int i = 0; i < listChromosomeFilePaths.length; i++) {
                    System.out.println("Processing Chromosome " + (i + 1));
                    String chromPath1 = listChromosomeFilePaths[i];
                    double[][] v1 = loadChromosomeFromFile(chromPath1);

                    for (int j = i + 1; j < listChromosomeFilePaths.length; j++) {
                        System.out.println("  Comparing with Chromosome " + (j + 1));
                        double[][] v2 = loadChromosomeFromFile(listChromosomeFilePaths[j]);

                        // ✅ 这里是唯一“算法层面”的改变：把 k 循环并行化
                        double freq = probabilityParallel(v1, v2, numberOfPoints, threads);

                        System.out.println("  Final link probability: " + freq);
                        out.println(freq * 100.0); // 与 C2 完全一致：输出百分比
                    }
                }
            }

            System.out.println("✅ Done: " + outPath.getAbsolutePath());
        }

        System.out.println("\n🎉 All folders finished.");
    }

    // =====================================================================
    // ✅ 与 C2 逻辑等价的概率计算：只是把 Monte Carlo(k) 并行化
    // 仍然调用 sampleTwoPoints + addedVector + compute_link_state
    // =====================================================================
    static double probabilityParallel(double[][] v1, double[][] v2, int numberOfPoints, int threads) {

        final AtomicInteger countLinked = new AtomicInteger(0);
        ExecutorService pool = Executors.newFixedThreadPool(threads);

        try {
            int chunk = (numberOfPoints + threads - 1) / threads;
            List<Future<?>> futs = new ArrayList<>();

            for (int tid = 0; tid < threads; tid++) {
                final int start = tid * chunk;
                final int end = Math.min(numberOfPoints, start + chunk);
                if (start >= end) break;

                futs.add(pool.submit(() -> {
                    int local = 0;
                    for (int k = start; k < end; k++) {

                        // ======= 以下三步保持 C2 完全一致 =======
                        double[][] randomPoints = sampleTwoPoints();
                        double[][] v1Aug = addedVector(v1, randomPoints[0]);
                        double[][] v2Aug = addedVector(v2, randomPoints[1]);
                        double linking = compute_link_state(v1Aug, v2Aug); // 0 or 1
                        // ======================================

                        if (linking != 0.0) local++;
                    }
                    countLinked.addAndGet(local);
                }));
            }

            for (Future<?> f : futs) {
                f.get();
            }

        } catch (Exception e) {
            throw new RuntimeException(e);
        } finally {
            pool.shutdown();
        }

        return ((double) countLinked.get()) / numberOfPoints;
    }

    // =====================================================================
    // ================== 以下：全部沿用 C2 原算法实现 ======================
    // =====================================================================

    public static double compute_link_state(double pathVect1[][], double pathVect2[][]) {

        double zERO = Math.pow(10, -8);

        double writhe_result = 0.0;
        double omegastar = 0.0;

        double r12[] = new double[3];
        double r13[] = new double[3];
        double r14[] = new double[3];
        double r23[] = new double[3];
        double r24[] = new double[3];
        double r34[] = new double[3];

        double n1[] = new double[3];
        double n2[] = new double[3];
        double n3[] = new double[3];
        double n4[] = new double[3];
        double n5[] = new double[3];

        double m1, m2, m3, m4;

        int i, j, k;

        double path[][] = pathVect1;
        int num_points = path.length;

        double path2[][] = pathVect2;
        int num_points2 = path2.length;

        for (i = 0; i < num_points; i++) {
            for (j = 0; j < num_points2; j++) {

                double pt1[] = new double[3];
                double pt2[] = new double[3];
                double pt3[] = new double[3];
                double pt4[] = new double[3];

                for (k = 0; k < 3; k++) {
                    pt1[k] = path[i][k];
                    pt2[k] = path[(i + 1) % num_points][k];

                    pt3[k] = path2[j][k];
                    pt4[k] = path2[(j + 1) % num_points2][k];

                    r12[k] = pt2[k] - pt1[k];
                    r13[k] = pt3[k] - pt1[k];
                    r14[k] = pt4[k] - pt1[k];
                    r23[k] = pt3[k] - pt2[k];
                    r24[k] = pt4[k] - pt2[k];
                    r34[k] = pt4[k] - pt3[k];
                }

                for (k = 0; k < 3; k++) {
                    n1[k] = Math.abs(vectorProduct(r13, r14, k)) <= zERO ? 0.0 : vectorProduct(r13, r14, k);
                    n2[k] = Math.abs(vectorProduct(r14, r24, k)) <= zERO ? 0.0 : vectorProduct(r14, r24, k);
                    n3[k] = Math.abs(vectorProduct(r24, r23, k)) <= zERO ? 0.0 : vectorProduct(r24, r23, k);
                    n4[k] = Math.abs(vectorProduct(r23, r13, k)) <= zERO ? 0.0 : vectorProduct(r23, r13, k);
                    n5[k] = Math.abs(vectorProduct(r34, r12, k)) <= zERO ? 0.0 : vectorProduct(r34, r12, k);
                }

                m1 = magnitude(n1);
                m2 = magnitude(n2);
                m3 = magnitude(n3);
                m4 = magnitude(n4);

                for (k = 0; k < 3; k++) {
                    if (m1 != 0) n1[k] /= m1;
                    if (m2 != 0) n2[k] /= m2;
                    if (m3 != 0) n3[k] /= m3;
                    if (m4 != 0) n4[k] /= m4;
                }

                omegastar = Math.asin(scalarProduct(n1, n2))
                        + Math.asin(scalarProduct(n2, n3))
                        + Math.asin(scalarProduct(n3, n4))
                        + Math.asin(scalarProduct(n4, n1));

                writhe_result += (omegastar * signum(scalarProduct(n5, r13))) / (4 * Math.PI);
            }
        }

        writhe_result = writhe_result * 2;

        if (Math.abs(writhe_result) < 0.1) {
            return 0.0;
        } else {
            return 1.0;
        }
    }

    // vector product component
    public static double vectorProduct(double vector1[], double vector2[], int i) {
        return (vector1[(i + 1) % 3] * vector2[(i + 2) % 3]
              - vector1[(i + 2) % 3] * vector2[(i + 1) % 3]);
    }

    // scalar product with clamp & near-zero as in C2
    public static double scalarProduct(double vector1[], double vector2[]) {

        double sc;
        double zERO = Math.pow(10, -8);

        sc = vector1[0] * vector2[0] +
             vector1[1] * vector2[1] +
             vector1[2] * vector2[2];

        if (sc > 1.0) return 1.0;
        else if (sc < -1.0) return -1.0;
        else if (Math.abs(sc) <= zERO) return 0;
        else return sc;
    }

    public static double magnitude(double vector[]) {
        return Math.sqrt(vector[0] * vector[0] +
                         vector[1] * vector[1] +
                         vector[2] * vector[2]);
    }

    public static int signum(double x) {
        double zERO = Math.pow(10, -8);
        if (x < -zERO) return -1;
        else if (x > zERO) return +1;
        else return 0;
    }

    public static void printVector(double vector[]) {
        System.out.printf("%f %f %f\n", vector[0], vector[1], vector[2]);
    }

    // I/O: whitespace-separated (same as C2)
    public static double[][] loadChromosomeFromFile(String fileToParse) throws FileNotFoundException, IOException {
        InputStream ips = new FileInputStream(fileToParse);
        InputStreamReader ipsr = new InputStreamReader(ips);
        ArrayList<double[]> tmpList = new ArrayList<>();

        try {
            BufferedReader br = new BufferedReader(ipsr);
            String ligne;

            while ((ligne = br.readLine()) != null) {
                String[] currentLineContent = ligne.trim().split("[,\\s]+");
                if (currentLineContent.length < 3) continue;

                double[] vector = new double[3];
                for (int i = 0; i < 3; i++) {
                    if (!currentLineContent[i].isEmpty()) {
                        double db = Double.valueOf(currentLineContent[i]);
                        vector[i] = db;
                    }
                }
                tmpList.add(vector);
            }

            ipsr.close();
            ips.close();

        } catch (IOException e) {
            System.out.println("Problem when reading data set");
        }

        double[][] matrice = new double[tmpList.size()][3];
        for (int i = 0; i < tmpList.size(); i++) {
            matrice[i] = tmpList.get(i);
        }
        return matrice;
    }

    // add a vector to a matrix (same as C2)
    public static double[][] addedVector(double[][] matrix, double[] vector) {
        double[][] mat = new double[matrix.length + 1][3];
        mat[0] = vector;
        for (int i = 0; i < matrix.length; i++) {
            mat[i + 1] = matrix[i];
        }
        return mat;
    }

    // sample two points on sphere (same as C2: uses Math.random)
    public static double[][] sampleTwoPoints() {
        double raduisOfSphere = 10000;
        double[][] points = new double[2][3];

        for (int i = 0; i < 2; ++i) {
            double z = (2.0 * Math.random() - 1.0) * raduisOfSphere;
            double t = 2.0 * Math.PI * Math.random();
            double r = Math.sqrt(raduisOfSphere * raduisOfSphere - z * z);

            points[i][0] = r * Math.cos(t);
            points[i][1] = r * Math.sin(t);
            points[i][2] = z;
        }
        return points;
    }

    // original probability() kept for reference (not used in parallel path)
    public static double probability(double vector[]) {
        double count = 0.0;
        for (int i = 0; i < vector.length; i++) {
            if (vector[i] != 0.0) count += 1.0;
        }
        return count / vector.length;
    }

    public static double[] centerOfMass(double[][] vector) {
        double sumXcoord = 0.0, sumYcoord = 0.0, sumZcoord = 0.0;
        double centerMass[] = new double[3];
        for (int i = 0; i < vector.length; i++) {
            sumXcoord += vector[i][0];
            sumYcoord += vector[i][1];
            sumZcoord += vector[i][2];
        }
        centerMass[0] = sumXcoord / vector.length;
        centerMass[1] = sumYcoord / vector.length;
        centerMass[2] = sumZcoord / vector.length;
        return centerMass;
    }

    public static double lagerRadius(double vector[], double matrix[][]) {
        double radius = 0.0;
        double maxRadius = 0.0;

        for (int i = 0; i < matrix.length; i++) {
            radius = 0.0;
            for (int j = 0; j < matrix[i].length; j++) {
                radius += (matrix[i][j] - vector[j]) * (matrix[i][j] - vector[j]);
            }
            radius = Math.sqrt(radius);

            if (maxRadius < radius) {
                maxRadius = radius;
            }
        }
        return maxRadius;
    }
}