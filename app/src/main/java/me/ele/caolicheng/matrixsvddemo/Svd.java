package me.ele.caolicheng.matrixsvddemo;

import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.SingularValueDecomposition;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by richardcao on 15/5/12.
 * svd构造和计算
 */
public class Svd {

    /**
     * 获取svd相关计算结果的方法：
     * list.get(0)是svd值的double[]，list中剩下的double[]是矩阵M*M的特征向量(V的列集合)
     *
     * @param data
     * @return
     */
    public static List<double[]> getSvd(double data[][]) throws Exception {
        DenseMatrix64F matA = buildMatrix(data);
        DebugLog.i(matA.toString());
        double[] result;
        ArrayList<double[]> resultList = new ArrayList<double[]>();
        SingularValueDecomposition<DenseMatrix64F> svd =
                DecompositionFactory.svd(matA.numRows, matA.numCols, false, true, true);
        if (!DecompositionFactory.decomposeSafe(svd, matA)) {
            throw new Exception("Decomposition failed");
        }
        DenseMatrix64F S = svd.getW(null);
        DenseMatrix64F V = svd.getV(null, false);
        result = new double[S.numRows];
        for (int i = 0; i < S.numRows; i++) {
            result[i] = formatDouble(S.get(i, i));
        }
        resultList.add(result);
        for (int i = 0; i < V.numCols; i++) {
            result = new double[V.numRows];
            for (int j = 0; j < V.numRows; j++) {
                result[j] = formatDouble(V.get(j, i));
            }
            resultList.add(result);
        }
        return resultList;
    }


    private static DenseMatrix64F buildMatrix(double data[][]) {
        return new DenseMatrix64F(data);
    }

    private static double formatDouble(double d) {
        DecimalFormat df = new DecimalFormat("#.####");
        return Double.parseDouble(df.format(d));
    }

    /**
     * 相似性度量计算
     *
     * @param svdA
     * @param svdB
     * @return
     */
    public static double getSimMeasure(List<double[]> svdA, List<double[]> svdB) throws Exception {
        double sum = 0;
        for (int i = 0; i < svdA.get(0).length; i++) {
            for (int j = 0; j < svdB.get(0).length; j++) {
                sum = sum + getWeight(i, svdA.get(0)) * getWeight(j, svdB.get(0))
                        * getDotProduct(svdA.get(i + 1), svdB.get(j + 1));
            }
        }
        return formatDouble(sum);
    }

    private static double getWeight(int i, double[] svd) {
        double sum = 0;
        for (double aSvd : svd) {
            sum = sum + Math.pow(aSvd, 2);
        }
        return (Math.pow(svd[i], 2)) / sum;
    }

    private static double getDotProduct(double[] a, double[] b) throws Exception {
        double sum = 0;
        if (a.length == b.length) {
            for (int i = 0; i < a.length; i++) {
                sum = sum + a[i] * b[i];
            }
        } else {
            throw new Exception("计算向量点积错误！");
        }
        return sum;
    }
}
