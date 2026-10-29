/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 *///============================================================================
// Name        : ComputeWrithe : to compute the linkingFrequency
// Author      : Maxime Guiffo Pouokam
// Version     :
// Copyright   : Copyrights by Owner
// Description : Java, Ansi-style
//============================================================================
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

/**
 *
 * @author maximepouokam
 */
public class ComputeWriteC2 {

 /**
     * @param args the command line arguments
     */
    
    static int numberOfSample;                 // 不给初值
    static String[] listChromosomeFilePaths;   // 不调用 create_Strings
    static final String RESULT_DIR = "C:/Users/weilu/Desktop/research william luo/HindIII+Mspl before FDR/java/RESULT";

// Create a string of paths  
    // === 把原来的 create_Strings(...) 替换成带 baseDir 的版本 ===
  public static String[] create_Strings(int simpleSize, String baseDir){
      String[] solution = new String[simpleSize];
      for (int i = 0; i < simpleSize; i++){
          int j = i + 1;
          // 用正斜杠，避免转义烦恼
          solution[i] = baseDir + "/" + j + ".txt";
      }
      return solution;
  }

  // === 新增：统计该目录下 1.txt, 2.txt, ... 连续存在到哪里 ===
  public static int countSequentialTxtFiles(String baseDir){
      int count = 0;
      while (true){
          java.io.File f = new java.io.File(baseDir + "/" + (count + 1) + ".txt");
          if (f.exists() && f.isFile()){
              count++;
          } else {
              break;
          }
      }
      return count;
  }

  // === 把原来的 main(...) 全部替换为这一版 ===
  public static void main(String[] args) throws FileNotFoundException, IOException {

      // 你要跑的四个输入目录（用 / 写法更稳）
      String[] baseDirs = new String[]{
          "C:\\Users\\weilu\\Downloads\\Panels (1)\\Panels\\Panel_1",
          "C:\\Users\\weilu\\Downloads\\Panels (1)\\Panels\\Panel_2"
          // "C:/Users/weilu/Desktop/research william luo/HindIII+Mspl before FDR/java/output_chroms_msplNMDS",
          // "C:/Users/weilu/Desktop/research william luo/HindIII+Mspl before FDR/java/output_chroms_msplPM1",
          // "C:/Users/weilu/Desktop/research william luo/HindIII+Mspl before FDR/java/output_chroms_msplPM2"
      };

      int numberOfPoints = 1000; // Monte Carlo 次数（可调）

      for (String baseDir : baseDirs){

          // 1) 统计该目录下有多少个按序编号的 txt（1.txt, 2.txt, ...）
          numberOfSample = countSequentialTxtFiles(baseDir);
          if (numberOfSample == 0){
              System.out.println("⚠️ 目录没有发现 1.txt, 2.txt... : " + baseDir);
              continue;
          }

          // 2) 生成该目录下的文件列表
          listChromosomeFilePaths = create_Strings(numberOfSample, baseDir);

          // 3) 为该目录准备一个独立的输出文件名（用目录最后一段做标签）
          String label = baseDir.replace('\\','/'); // 统一分隔符
          int lastSlash = label.lastIndexOf('/');
          label = (lastSlash >= 0) ? label.substring(lastSlash + 1) : label; // 取最后一段

          // 统一结果目录
          java.io.File resultDir = new java.io.File(RESULT_DIR);
          resultDir.mkdirs();


          java.io.File outPath = new java.io.File(resultDir, label + "_linking_frequency.txt");

          System.out.println("\n========================");
          System.out.println("▶ Running on folder: " + baseDir);
          System.out.println("▶ Samples detected:  " + numberOfSample);
          System.out.println("▶ Output file:       " + outPath.getAbsolutePath());
          System.out.println("========================");

          // 开始写文件到 RESULT 目录
          try (FileWriter fw = new FileWriter(outPath);
              PrintWriter out = new PrintWriter(new BufferedWriter(fw))) {

              for (int i = 0; i < listChromosomeFilePaths.length; i++){
                  System.out.println("Processing Chromosome " + (i+1));
                  String chromPath1 = listChromosomeFilePaths[i];
                  double[][] v1 = loadChromosomeFromFile(chromPath1);

                  for (int j = i + 1; j < listChromosomeFilePaths.length; j++){
                      System.out.println("  Comparing with Chromosome " + (j+1));
                      double[][] v2 = loadChromosomeFromFile(listChromosomeFilePaths[j]);

                      double[] vectorOfBinary = new double[numberOfPoints];

                      for (int k = 0; k < numberOfPoints; k++){
                          double[][] randomPoints = sampleTwoPoints();
                          double[][] v1Aug = addedVector(v1, randomPoints[0]);
                          double[][] v2Aug = addedVector(v2, randomPoints[1]);

                          double linking = compute_link_state(v1Aug, v2Aug); // 0 or 1
                          vectorOfBinary[k] = linking;
                      }

                      double freq = probability(vectorOfBinary); // 0~1
                      System.out.println("  Final link probability: " + freq);
                      out.println(freq * 100.0); // 与你原始逻辑保持一致：输出百分比
                  }
              }
          }

          System.out.println("✅ Done: " + outPath.getAbsolutePath());
      }

      System.out.println("\n🎉 All folders finished.");
  }

    //// Calculates the writhe using the Gauss integral
// described in Klein and Langowski 2000 [KL]
// also calculates ACN by using |omega*|

    

public static double compute_link_state(double pathVect1[][], double pathVect2[][])
{
    
 // Long start = System.currentTimeMillis() ;    // get time 
  double zERO=Math.pow(10, -8);;

  double writhe_result = 0.0; // answers
 // double acn_result = 0;
  double omegastar = 0.0; // omega star as in KL
  double term; // term in the sum (omega_ij)
  // disp. vectors
  double r12[]=new double[3];
  double r13[]=new double[3];
  double r14[]=new double[3];
  double r23[]=new double[3];
  double r24[]=new double[3];
  double r34[]=new double[3];
          
   // n vectors used to calculate omega_ij
  double n1[]=new double[3];
  double n2[]=new double[3];
  double n3[]=new double[3];
  double n4[]=new double[3];
                              
  double n5[]=new double[3];// used to calculate the sign (n5 = r34 x r12)
  double m1, m2, m3, m4; // magnitudes of the cross products

  int i,j,k; // counters
  
  double path[][]= pathVect1;
  //{{0,0,0},{2 , 0 ,0 },{2 , 0, 4},{0, 0 , 4}};
  //{{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,0}};
  //{{0,0,0},{1,0,0},{1,1,0},{1,1,1},{1,0,1},{0,0,1},{0,0,0}}; //path of the first curve 
  
 // {{1,1,1},{2,1,1},{2,2,1},{2,2,2},{2,1,2},{1,1,2},{1,1,1}};
  int num_points=path.length; 
 //double path2[][]={{0.75,0.75,0.75},
 // {1.75,0.75,0.75},{1.75,1.75,0.75},{1.75,1.75,1.75},{1.75,0.75,1.75},{0.75,0.75,1.75},};// path of the second curve
 double  path2[][]=pathVect2;
         //{{0.75,0.75,0.75},{1.75,0.75,0.75},{1.75,1.75,0.75},{0.75,1.75,0.75},{0.75,0.75,0.75}};
  
 //{{1,-1,2}, {3,-1,2}, {3,3,2},{1,3,2}};


  int num_points2=path2.length ;
  

  for(i=0;i<num_points;i++) // first sum (note off by 1 from K&L)
    {
    for(j=0;j<num_points2;j++) // second sum
      {

      // Calculate the displacement vectors.
      // segment i runs from pt 1 to pt 2.  segment j from pt 3 to pt 4
     
      double pt1[]=new double[3];
      double pt2[]=new double[3];
      double pt3[]=new double[3];
      double pt4[]=new double[3];
      for(k=0;k<3;k++)
        {
        pt1[k] = path[i][k];
        pt2[k] = path[(i+1)%num_points][k];

        pt3[k] = path2[j][k];
        pt4[k] = path2[(j+1)%num_points2][k];

        r12[k] = pt2[k] - pt1[k];
        r13[k] = pt3[k] - pt1[k];
        r14[k] = pt4[k] - pt1[k];
        r23[k] = pt3[k] - pt2[k];
        r24[k] = pt4[k] - pt2[k];
        r34[k] = pt4[k] - pt3[k];
        }

      // calculate the unnormalized cross products and fix round-off error
      for(k=0;k<3;k++)
        {
        n1[k] = Math.abs(vectorProduct(r13,r14,k)) <= zERO ? 0.0 :
                     vectorProduct(r13,r14,k);
        n2[k] = Math.abs(vectorProduct(r14,r24,k)) <= zERO ? 0.0 :
                     vectorProduct(r14,r24,k);
        n3[k] = Math.abs(vectorProduct(r24,r23,k)) <= zERO ? 0.0 :
                     vectorProduct(r24,r23,k);
        n4[k] = Math.abs(vectorProduct(r23,r13,k)) <= zERO ? 0.0 :
                     vectorProduct(r23,r13,k);
        //n5 used in calculating the sign of omegastar
        n5[k] = Math.abs(vectorProduct(r34,r12,k)) <= zERO ? 0.0 :
                     vectorProduct(r34,r12,k);
        }

      // normalize these vectors
      m1 = magnitude(n1);
      m2 = magnitude(n2);
      m3 = magnitude(n3);
      m4 = magnitude(n4);
      for(k=0;k<3;k++)
        {
        if (m1!=0) {n1[k] /= m1;}
        
        if (m2!=0) {n2[k] /= m2;}
        if (m3!=0) {n3[k] /= m3;}
        if (m4!=0) {n4[k] /= m4;}
        }

      // Calculate omega* since n are normalized, can get angles
      // directly from the dot product.
      omegastar=Math.asin(scalarProduct(n1,n2)) +
                Math.asin(scalarProduct(n2,n3)) +
                Math.asin(scalarProduct(n3,n4)) +
                Math.asin(scalarProduct(n4,n1));

      // add the sign and update sum
      writhe_result += (omegastar*signum(scalarProduct(n5,r13)))/(4*Math.PI);
    //  acn_result += omegastar/(4*Math.PI); // always positive

      }
    }

  writhe_result = writhe_result*2;  //Multiply final sum by 2
 // acn_result = acn_result*2;

 // return;
  
 // Long time = System.currentTimeMillis() - start ;   //get time
 // System.out.println("exec time = " + time);         // get time
  if (Math.abs(writhe_result) <0.1){
      //writhe_result=0;
      return 0.0;
  }   
       
  else {
      return 1.0; }
  
  //System.out.println(writhe_result);
 // System.out.println(acn_result);
 //return writhe_result;
}
// ----------------------------------------------------------------- COMPUTE LINKING NUMBER 
    

public static double compute_linking(double pathVect1[][], double pathVect2[][])
{
    
 // Long start = System.currentTimeMillis() ;    // get time 
  double zERO=Math.pow(10, -8);;

  double writhe_result = 0.0; // answers
  //double acn_result = 0.0;
  double omegastar = 0.0; // omega star as in KL
  double term; // term in the sum (omega_ij)
  // disp. vectors
  double r12[]=new double[3];
  double r13[]=new double[3];
  double r14[]=new double[3];
  double r23[]=new double[3];
  double r24[]=new double[3];
  double r34[]=new double[3];
          
   // n vectors used to calculate omega_ij
  double n1[]=new double[3];
  double n2[]=new double[3];
  double n3[]=new double[3];
  double n4[]=new double[3];
                              
  double n5[]=new double[3];// used to calculate the sign (n5 = r34 x r12)
  double m1, m2, m3, m4; // magnitudes of the cross products

  int i,j,k; // counters
  
  double path[][]= pathVect1;
  //{{0,0,0},{2 , 0 ,0 },{2 , 0, 4},{0, 0 , 4}};
  //{{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,0}};
  //{{0,0,0},{1,0,0},{1,1,0},{1,1,1},{1,0,1},{0,0,1},{0,0,0}}; //path of the first curve 
  
 // {{1,1,1},{2,1,1},{2,2,1},{2,2,2},{2,1,2},{1,1,2},{1,1,1}};
  int num_points=path.length; 
 //double path2[][]={{0.75,0.75,0.75},
 // {1.75,0.75,0.75},{1.75,1.75,0.75},{1.75,1.75,1.75},{1.75,0.75,1.75},{0.75,0.75,1.75},};// path of the second curve
 double  path2[][]=pathVect2;
         //{{0.75,0.75,0.75},{1.75,0.75,0.75},{1.75,1.75,0.75},{0.75,1.75,0.75},{0.75,0.75,0.75}};
  
 //{{1,-1,2}, {3,-1,2}, {3,3,2},{1,3,2}};


  int num_points2=path2.length ;
  
//  System.out.println(num_points2);
  //System.out.println(num_points);
  
  

  for(i=0;i<num_points;i++) // first sum (note off by 1 from K&L)
    {
    for(j=0;j<num_points2;j++) // second sum
      {

      // Calculate the displacement vectors.
      // segment i runs from pt 1 to pt 2.  segment j from pt 3 to pt 4
     
      double pt1[]=new double[3];
      double pt2[]=new double[3];
      double pt3[]=new double[3];
      double pt4[]=new double[3];
      for(k=0;k<3;k++)
        {
        pt1[k] = path[i][k];
        pt2[k] = path[(i+1)%num_points][k];

        pt3[k] = path2[j][k];
        pt4[k] = path2[(j+1)%num_points2][k];

        r12[k] = pt2[k] - pt1[k];
        r13[k] = pt3[k] - pt1[k];
        r14[k] = pt4[k] - pt1[k];
        r23[k] = pt3[k] - pt2[k];
        r24[k] = pt4[k] - pt2[k];
        r34[k] = pt4[k] - pt3[k];
        }

      // calculate the unnormalized cross products and fix round-off error
      for(k=0;k<3;k++)
        {
        n1[k] = Math.abs(vectorProduct(r13,r14,k)) <= zERO ? 0.0 :
                     vectorProduct(r13,r14,k);
        n2[k] = Math.abs(vectorProduct(r14,r24,k)) <= zERO ? 0.0 :
                     vectorProduct(r14,r24,k);
        n3[k] = Math.abs(vectorProduct(r24,r23,k)) <= zERO ? 0.0 :
                     vectorProduct(r24,r23,k);
        n4[k] = Math.abs(vectorProduct(r23,r13,k)) <= zERO ? 0.0 :
                     vectorProduct(r23,r13,k);
        //n5 used in calculating the sign of omegastar
        n5[k] = Math.abs(vectorProduct(r34,r12,k)) <= zERO ? 0.0 :
                     vectorProduct(r34,r12,k);
        }

      // normalize these vectors
      m1 = magnitude(n1);
      m2 = magnitude(n2);
      m3 = magnitude(n3);
      m4 = magnitude(n4);
      for(k=0;k<3;k++)
        {
        if (m1!=0) {n1[k] /= m1;}
        
        if (m2!=0) {n2[k] /= m2;}
        if (m3!=0) {n3[k] /= m3;}
        if (m4!=0) {n4[k] /= m4;}
        }

      // Calculate omega* since n are normalized, can get angles
      // directly from the dot product.
      omegastar=Math.asin(scalarProduct(n1,n2)) +
                Math.asin(scalarProduct(n2,n3)) +
                Math.asin(scalarProduct(n3,n4)) +
                Math.asin(scalarProduct(n4,n1));

      // add the sign and update sum
      writhe_result += (omegastar*signum(scalarProduct(n5,r13)))/(4*Math.PI);
     // acn_result += omegastar/(4*Math.PI); // always positive

      }
    }

  //writhe_result = writhe_result;  //Multiply final sum by 2
  
  //writhe_result = writhe_result;  //Multiply final sum by 2
 // acn_result = acn_result*2;

 // return;
  
 // Long time = System.currentTimeMillis() - start ;   //get time
 // System.out.println("exec time = " + time);         // get time
  //if (Math.abs(writhe_result) <0.1){
      //writhe_result=0;
    //  return 0.0;
  //}   
       
  //else {
    //  return 1.0; }
  
  //System.out.println(writhe_result);
// System.out.println(acn_result);//
 return writhe_result;
}



//----------------------------------------------------------------------------------------- END LINKING NUMBER 

// ---------------------------------- LINKING OPEN 

public static double compute_linking_open(double pathVect1[][], double pathVect2[][])
{
    
 // Long start = System.currentTimeMillis() ;    // get time 
  double zERO=Math.pow(10, -8);;

  double writhe_result = 0.0; // answers
 // double acn_result = 0;
  double omegastar = 0.0; // omega star as in KL
  double term; // term in the sum (omega_ij)
  // disp. vectors
  double r12[]=new double[3];
  double r13[]=new double[3];
  double r14[]=new double[3];
  double r23[]=new double[3];
  double r24[]=new double[3];
  double r34[]=new double[3];
          
   // n vectors used to calculate omega_ij
  double n1[]=new double[3];
  double n2[]=new double[3];
  double n3[]=new double[3];
  double n4[]=new double[3];
                              
  double n5[]=new double[3];// used to calculate the sign (n5 = r34 x r12)
  double m1, m2, m3, m4; // magnitudes of the cross products

  int i,j,k; // counters
  
  double path[][]= pathVect1;
  //{{0,0,0},{2 , 0 ,0 },{2 , 0, 4},{0, 0 , 4}};
  //{{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,0}};
  //{{0,0,0},{1,0,0},{1,1,0},{1,1,1},{1,0,1},{0,0,1},{0,0,0}}; //path of the first curve 
  
 // {{1,1,1},{2,1,1},{2,2,1},{2,2,2},{2,1,2},{1,1,2},{1,1,1}};
  int num_points=path.length  ; 
 //double path2[][]={{0.75,0.75,0.75},
 // {1.75,0.75,0.75},{1.75,1.75,0.75},{1.75,1.75,1.75},{1.75,0.75,1.75},{0.75,0.75,1.75},};// path of the second curve
 double  path2[][]=pathVect2;
         //{{0.75,0.75,0.75},{1.75,0.75,0.75},{1.75,1.75,0.75},{0.75,1.75,0.75},{0.75,0.75,0.75}};
  
 //{{1,-1,2}, {3,-1,2}, {3,3,2},{1,3,2}};


  int num_points2=path2.length ;
  

  for(i=0;i<num_points -1;i++) // first sum (note off by 1 from K&L)
    {
    for(j=0;j<num_points2-1;j++) // second sum
      {

      // Calculate the displacement vectors.
      // segment i runs from pt 1 to pt 2.  segment j from pt 3 to pt 4
     
      double pt1[]=new double[3];
      double pt2[]=new double[3];
      double pt3[]=new double[3];
      double pt4[]=new double[3];
      
      for(k=0;k<3;k++)
        {
        pt1[k] = path[i][k];
        pt2[k] = path[(i+1)%num_points][k];

        pt3[k] = path2[j][k];
        pt4[k] = path2[(j+1)%num_points2][k];

        r12[k] = pt2[k] - pt1[k];
        r13[k] = pt3[k] - pt1[k];
        r14[k] = pt4[k] - pt1[k];
        r23[k] = pt3[k] - pt2[k];
        r24[k] = pt4[k] - pt2[k];
        r34[k] = pt4[k] - pt3[k];
        }

      // calculate the unnormalized cross products and fix round-off error
      for(k=0;k<3;k++)
        {
        n1[k] = Math.abs(vectorProduct(r13,r14,k)) <= zERO ? 0.0 :
                     vectorProduct(r13,r14,k);
        n2[k] = Math.abs(vectorProduct(r14,r24,k)) <= zERO ? 0.0 :
                     vectorProduct(r14,r24,k);
        n3[k] = Math.abs(vectorProduct(r24,r23,k)) <= zERO ? 0.0 :
                     vectorProduct(r24,r23,k);
        n4[k] = Math.abs(vectorProduct(r23,r13,k)) <= zERO ? 0.0 :
                     vectorProduct(r23,r13,k);
        //n5 used in calculating the sign of omegastar
        n5[k] = Math.abs(vectorProduct(r34,r12,k)) <= zERO ? 0.0 :
                     vectorProduct(r34,r12,k);
        }

      // normalize these vectors
      m1 = magnitude(n1);
      m2 = magnitude(n2);
      m3 = magnitude(n3);
      m4 = magnitude(n4);
      
      for(k=0;k<3;k++)
        {
        if (m1!=0) {n1[k] /= m1;}
        
        if (m2!=0) {n2[k] /= m2;}
        if (m3!=0) {n3[k] /= m3;}
        if (m4!=0) {n4[k] /= m4;}
        }

      // Calculate omega* since n are normalized, can get angles
      // directly from the dot product.
      omegastar=Math.asin(scalarProduct(n1,n2)) +
                Math.asin(scalarProduct(n2,n3)) +
                Math.asin(scalarProduct(n3,n4)) +
                Math.asin(scalarProduct(n4,n1));

      // add the sign and update sum
      writhe_result += (omegastar*signum(scalarProduct(n5,r13)))/(4*Math.PI);
     // System.out.println(writhe_result);
    //  acn_result += omegastar/(4*Math.PI); // always positive

      }
    }

  //writhe_result = writhe_result;  //Multiply final sum by 2
 // System.out.println(writhe_result);
 // acn_result = acn_result*2;

 // return;
  
 // Long time = System.currentTimeMillis() - start ;   //get time
 // System.out.println("exec time = " + time);         // get time
  //if (Math.abs(writhe_result) <0.1){
      //writhe_result=0;
   //   return 0.0;
  //}   
       
  //else {
  //    return 1.0; }
  
  //System.out.println(writhe_result);
 // System.out.println(acn_result);
 return writhe_result;
}



//-------------------------------------END LINKING OPEN 

    // Calculates the ith entry of the vector product
public static double vectorProduct(double vector1[], double vector2[], int i)
{

  return (vector1[(i+1)%3]*vector2[(i+2)%3] -
                    vector1[(i+2)%3]*vector2[(i+1)%3]);
}

// Calculates the Scalar product of two vectors
public static double scalarProduct(double vector1[], double vector2[])
{

  double sc;
  double zERO=Math.pow(10, -8);;

  sc = vector1[0]*vector2[0] +
       vector1[1]*vector2[1] +
       vector1[2]*vector2[2];

  // Must keep range of scalar product within normalized bounds [-1..1]
  if(sc > 1.0)
    {
    return(1.0);
    }
  else if(sc < -1.0)
    {
    return(-1.0);
    }
  else if(Math.abs(sc) <= zERO)
    {
    return(0);
    }
  else{
    return sc;
}
}
    // Calculates the magnitude of a vector
// This does not limit the bounds like ScalarProduct()
    
public static double magnitude(double vector[])// vector[3]
{
  return Math.sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]);
}

// signum
// this ZERO is used to determine if a double is zero or not (because of
// rounding errors, a number that is zero often ends up being slightly
// nonzero)
public static int signum(double x)
        
{
    double zERO=Math.pow(10, -8);;
  if(x < -zERO)
    {
    return -1;
    }
  else if (x > zERO)
    {
    return +1;
    }
  else
  {
    return 0;
  }
}

    // Prints the vector
public static void printVector(double vector[])// showed vector[3]
{
  System.out.printf("%f %f %f\n", vector[0], vector[1], vector[2]);
}


 public static double[][] loadChromosomeFromFile(String fileToParse) throws FileNotFoundException, IOException {
        InputStream ips = new FileInputStream(fileToParse);
        InputStreamReader ipsr = new InputStreamReader(ips);
        ArrayList<double[]> tmpList = new ArrayList();
        try {
            BufferedReader br = new BufferedReader(ipsr);
            String ligne ;
            
            while ((ligne = br.readLine()) != null ) {

                // String[] currentLineContent =  ligne.split(",", -0); // this one for comma separated files
                String[] currentLineContent = ligne.trim().split("\\s+"); // this one for space separated files

                       // ligne.split("	", -1);
                        //ligne.split("", -1);
                                            // ligne.split("  ", -1);
               //  line.split("    ");
               
                // System.out.println( currentLineContent);
             

                double[] vector = new double[3];
                for (int i = 0; i < 3; i++) {
                    if (!currentLineContent[i].isEmpty()) {
                        
                        double db = Double.valueOf(currentLineContent[i]);
                        //System.out.println(db);           
                        vector[i] = db ;
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
            double[] xyz = tmpList.get(i);
            matrice[i] = xyz;
        }

//    for (int i = 0; i < matrice.length; i++) {
//          double[] xyz = matrice[i];
//         System.out.println(xyz[0] + "," + xyz[1] + "," + xyz[2]);
//     }
        return matrice;
    }
 
 // added a vector to a matrix
 public static double[][] addedVector(double[][] matrix, double[] vector){
     
     double[][] mat= new double[matrix.length+1][3];
     
     mat[0]=vector;
     
     
     for (int i=0; i< matrix.length;i++){
         mat[i+1]=matrix[i];
         
         
     }
     
     return mat;
  }
 
 // compute the probability  non zero numbers in a given vector
  public static double probability( double vector[]){
           
           double count=0.0;
           
           for (int i =0; i< vector.length ; i++){
               
               if (vector[i] !=0.0 ){
                   count+= 1.0;
               }
              }
           return count/vector.length;
       }
  // compute the center of mass of a given Nx3 vector
  
       public static double[] centerOfMass(double[][] vector){
           
           double sumXcoord=0.0;
           double sumYcoord=0.0;
           double sumZcoord=0.0;
           double centerMass[]= new double[3];
           for (int i=0; i<vector.length;i++){
                sumXcoord +=vector[i][0];
                sumYcoord +=vector[i][1];
                sumZcoord +=vector[i][2];
               
           }
           
           centerMass[0]=sumXcoord/vector.length;
           centerMass[1]=sumYcoord/vector.length;
           centerMass[2]=sumZcoord/vector.length;

            
           return centerMass;
           
       }
       
       // compute the lagerst radius of a sphere that emcompasses a set of points
       // given by a matrix and center at vector
       
       public static double lagerRadius(double vector[], double matrix[][]){
           
           double radius=0.0;
           double maxRadius=0.0;
           
            for (int i=0; i< matrix.length;i++ ){
               
               for (int j=0; j< matrix[i].length;j++){
                   
                   radius += (matrix[i][j]-vector[j])*(matrix[i][j]-vector[j]);
               }
               radius=Math.sqrt(radius);
               
               if (maxRadius < radius){
                   
                   maxRadius=radius;
                   
                   radius=0.0;
               }
               
           }
           return maxRadius;
       }
       
       
       
       
       public static double[][] sampleTwoPoints(){
           double raduisOfSphere =  10000 ; // 100000
           //double[] transLaVec= {108.3, 99, 99};
           
           double[][] points=new double[2][3];
 
	for(int i=0; i<2; ++i)
	{
		double z = (2.0*Math.random() - 1.0)*raduisOfSphere;
		double t = 2.0*Math.PI*Math.random();
		double r = Math.sqrt(raduisOfSphere*raduisOfSphere-z*z);
		
                
                points[i][0]=r*Math.cos(t);
                points[i][1]=r*Math.sin(t);
                points[i][2]=z;
                
	}
        
        //System.out.printf("%f %f %f \n", points[1][0],points[1][1],points[1][2]);
        return points;
}
       
       
       
	
}
