
using System;

namespace EikonalSolver;

public class Eikonal2D
{
    int m;
    int n;
    float[][] T;
    float[][] S;
    float dx;
    float dy;
    float dx2;
    float dy2;
    float dx2dy2;
    float dr2;
    int nBand = 0;
    int nBandMax = 0;
    int[][] band;
    int[] mx;
    int[] px;
    int[] my;
    int[] py;
    int xI;
    int xJ;
    int[] iBand;
    int[] jBand;
    int count;

    public Eikonal2D(float[][] T, float[][] S, float dx, float dy)
    {
        this.T = T;
        this.S = S;
        this.dx = dx;
        this.dy = dy;

        m = T.length;
        n = T[0].length;

        band = new int[m][n];
        mx = new int[m];
        px = new int[m];
        my = new int[n];
        py = new int[n];

        for (int i = 0; i < m; i++)
        {
            mx[i] = i - 1;
            px[i] = i + 1;
        }
        mx[0] = 1;
        px[m - 1] = m - 2;

        for (int i = 0; i < n; i++)
        {
            my[i] = i - 1;
            py[i] = i + 1;
        }
        my[0] = 1;
        py[n - 1] = n - 2;

        dx2 = dx * dx;
        dy2 = dy * dy;
        dx2dy2 = dx2 * dy2;
        dr2 = dx2 + dy2;
    }

    private void UpdateTravelTime(int i, int j)
    {
        float Sij = S[i][j];
        float a2;
        float a;
        float a2;
        if (T[mx[i]][j] < T[px[i]][j])
        {
            float a = T[mx[i]][j];
            a2 = T[mx[mx[i]]][j];
        }
        else
        {
            a = T[px[i]][j];
            a2 = T[px[px[i]]][j];
        }
        float b2;
        float b;
        float b2;
        if (T[i][my[j]] < T[i][py[j]])
        {
            float b = T[i][my[j]];
            b2 = T[i][my[my[j]]];
        }
        else
        {
            b = T[i][py[j]];
            b2 = T[i][py[py[j]]];
        }
        float ta = a + Sij * dx;
        float tb = b + Sij * dy;
        float tMin = ta < tb ? ta : tb;
        bool aOK = a < tMin;
        bool bOK = b < tMin;
        if (aOK && bOK)
        {
            float EE;
            float AA;
            float BB;
            float DD;
            float EE;
            if (a2 < a && b2 < b)
            {
                float AA = 4 * a - a2;
                float BB = 4 * b - b2;
                float DD = 3 * (dx2 + dy2);
                EE = dx2dy2 * (dr2 * Sij * Sij * 4 - (AA - BB) * (AA - BB));
            }
            else
            {
                AA = a;
                BB = b;
                DD = dx2 + dy2;
                EE = dx2dy2 * (DD * Sij * Sij - (AA - BB) * (AA - BB));
            }

            if (EE > 0)
                tMin = (AA * dy2 + BB * dx2 + MathF.Sqrt(EE)) / DD;
        }
        T[i][j] = tMin;
        UpdateBand(i, j);
    }

    private void UpdateBand(int i, int j)
    {
        if (band[i][j] == 0)
        {
            nBand += 1;
            if (nBand > nBandMax) ExtendHeap();
            iBand[nBand] = i;
            jBand[nBand] = j;
            band[i][j] = nBand;
        }
        UpHeap(band[i][j]);
    }

    private void MinTravelTimeInBand()
    {
        xI = iBand[1];
        xJ = jBand[1];
        band[xI][xJ] = -1;
        iBand[1] = iBand[nBand];
        jBand[1] = jBand[nBand];
        band[iBand[1]][jBand[1]] = 1;
        nBand -= 1;
        DownHeap(1);
    }

    public void UpHeap(int item)
    {
        int p = item / 2;
        int i = iBand[item];
        int j = jBand[item];
        int k = item;
        while (((k > 1 ? 1 : 0) & (T[iBand[p]][jBand[p]] > T[i][j] ? 1 : 0)) != 0)
        {
            iBand[k] = iBand[p];
            jBand[k] = jBand[p];
            band[iBand[k]][jBand[k]] = k;
            k = p;
            p = k / 2;
        }
        iBand[k] = i;
        jBand[k] = j;
        band[i][j] = k;
    }

    public void DownHeap(int item)
    {
        int i = iBand[item];
        int j = jBand[item];
        int k = item;
        while (k <= nBand / 2)
        {
            int m = k + k;
            int n = m + 1;
            if (m < nBand && T[iBand[m]][jBand[m]] > T[iBand[n]][jBand[n]]) m = n;

            if (T[i][j] < T[iBand[m]][jBand[m]]) break;

            iBand[k] = iBand[m];
            jBand[k] = jBand[m];
            band[iBand[k]][jBand[k]] = k;
            k = m;
        }
        iBand[k] = i;
        jBand[k] = j;
        band[i][j] = k;
    }

    public void ExtendHeap()
    {
        int[] iTemp = iBand;
        int[] jTemp = jBand;
        int oldMax = nBandMax;

        if (nBandMax == 0)
            nBandMax = n + m;
        else
            nBandMax = 2 * nBandMax;

        iBand = new int[nBandMax + 1];
        jBand = new int[nBandMax + 1];

        if (oldMax > 0)
        {
            Array.Copy(iTemp, 1, iBand, 1, oldMax);
            Array.Copy(jTemp, 1, jBand, 1, oldMax);
        }
    }

    public void Solve()
    {
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
                band[i][j] = float.IsNaN(T[i][j]) ? 0 : -1;
        }
        for (i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (band[i][j] == -1)
                {
                    if (band[mx[i]][j] != -1) UpdateTravelTime(mx[i], j);
                    if (band[px[i]][j] != -1) UpdateTravelTime(px[i], j);
                    if (band[i][my[j]] != -1) UpdateTravelTime(i, my[j]);
                    if (band[i][py[j]] != -1) UpdateTravelTime(i, py[j]);
                }
            }
        }

        while (nBand > 0)
        {
            count += 1;
            MinTravelTimeInBand();
            i = xI; 
            
            int j = xJ;
            if (band[mx[i]][j] != -1) UpdateTravelTime(mx[i], j);
            if (band[px[i]][j] != -1) UpdateTravelTime(px[i], j);
            if (band[i][my[j]] != -1) UpdateTravelTime(i, my[j]);
            if (band[i][py[j]] != -1) UpdateTravelTime(i, py[j]);
        }
    }
}