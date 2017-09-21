package org.panda.tfactivity;

import org.panda.resource.tcga.ExpressionReader;
import org.panda.utility.statistics.FDR;
import org.panda.utility.statistics.Summary;
import org.panda.utility.statistics.TTest;

import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Ozgun Babur
 */
public class TCGAExpressionProvider implements DiscreteExpressionProvider
{
	ExpressionReader er;
	TwoGroupsSampleNames two;
	double pvalThr;

	Map<String, Integer> cache;

	public TCGAExpressionProvider(String file) throws FileNotFoundException
	{
		er = new ExpressionReader(file);
		cache = new HashMap<>();
	}

	public void setTwo(TwoGroupsSampleNames two)
	{
		this.two = two;
		cache.clear();
	}

	public void setPvalThr(double pvalThr)
	{
		this.pvalThr = pvalThr;
		cache.clear();
	}

	public void setFDRThr(double thr)
	{
		setPvalThr(getPvalThrForGivenFDR(thr));
	}

	@Override
	public Integer getChange(String gene)
	{
		if (!cache.containsKey(gene))
		{
			Integer change = calcChange(gene);
			cache.put(gene, change);
		}
		return cache.get(gene);
	}

	private Integer calcChange(String gene)
	{
		double[] ctrl = er.getGeneAlterationArray(gene, two.getControl());
		double[] test = er.getGeneAlterationArray(gene, two.getTest());
		Double pval = getPval(ctrl, test);

		if (pval == null || pval > pvalThr) return 0;

		if (Summary.mean(test) > Summary.mean(ctrl)) return 1;
		return -1;
	}

	private Double getPval(String gene)
	{
		double[] ctrl = er.getGeneAlterationArray(gene, two.getControl());
		double[] test = er.getGeneAlterationArray(gene, two.getTest());

		return getPval(ctrl, test);
	}

	private Double getPval(double[] ctrl, double[] test)
	{
		if (ctrl == null || test == null) return null;

		double pval = TTest.getPValOfMeanDifference(ctrl, test);
		if (Double.isNaN(pval)) return null;
		return pval;
	}

	private double getPvalThrForGivenFDR(double fdr)
	{
		Map<String, Double> map = new HashMap<>();
		for (String gene : er.getGenes())
		{
			Double p = getPval(gene);
			if (p != null && !p.isNaN())
			{
				map.put(gene, p);
			}
		}

		double thr = FDR.getPValueThreshold(map, null, fdr);
		System.out.println("pval thr matching given fdr = " + thr);
		return thr;
	}
}
