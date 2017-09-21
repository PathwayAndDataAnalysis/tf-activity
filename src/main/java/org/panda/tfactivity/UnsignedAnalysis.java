package org.panda.tfactivity;

import org.panda.utility.Progress;
import org.panda.utility.statistics.FDR;
import org.panda.utility.statistics.FishersExactTest;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class UnsignedAnalysis
{
	/**
	 * Unsigned transcriptional network.
	 */
	Map<String, Set<String>> network;

	/**
	 * Expression data provider.
	 */
	DiscreteExpressionProvider expProv;

	public UnsignedAnalysis(Map<String, Set<String>> network, DiscreteExpressionProvider expProv)
	{
		this.network = network;
		this.expProv = expProv;
	}

	public void run(String outFile, double fdrThr) throws IOException
	{
		Map<String, Double> pValues = getActivityPValues(network, expProv);
		writeResults(outFile, network, expProv, pValues, fdrThr);
	}

	private Map<String, Double> getActivityPValues(Map<String, Set<String>> network,
		DiscreteExpressionProvider expProv)
	{
		Map<String, Double> pvals = new HashMap<>();

		int size = (int) network.values().stream().flatMap(Collection::stream).distinct().count();

		int featured = (int) network.values().stream().flatMap(Collection::stream)
			.filter(t -> expProv.getChange(t) != 0).distinct().count();

		for (String tf : network.keySet())
		{
			Set<String> targets = network.get(tf);

			int selected = targets.size();
			int featuredSelected = (int) targets.stream().filter(t -> expProv.getChange(t) != 0).count();

			double pval = FishersExactTest.calcEnrichmentPval(size, featured, selected, featuredSelected);
			pvals.put(tf, pval);
		}

		return pvals;
	}

	private Map<String, List<String>> getActivitySupporterGenes(Map<String, Set<String>> network,
		DiscreteExpressionProvider expProv)
	{
		Map<String, List<String>> evidenceMap = new HashMap<>();

		for (String tf : network.keySet())
		{
			List<String> list = new ArrayList<>();

			Set<String> set = network.get(tf);

			for (String target : set)
			{
				if (expProv.getChange(target) != 0)
				{
					list.add(target);
				}
			}

			Collections.sort(list);
			evidenceMap.put(tf, list);
		}

		return evidenceMap;
	}

	private void writeResults(String file, Map<String, Set<String>> network,
		DiscreteExpressionProvider expProv, Map<String, Double> pMap, double fdrThr) throws IOException
	{
		Map<String, Double> qMap = FDR.getQVals(pMap, null);

		List<String> select = FDR.select(pMap, null, fdrThr);

		Map<String, List<String>> support = getActivitySupporterGenes(network, expProv);

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(file));
		writer.write("TF\tP-val\tQ-val\tSupporting targets");

		for (String tf : select)
		{
			writer.write("\n" + tf + "\t" + pMap.get(tf) + "\t" + qMap.get(tf) + "\t" + support.get(tf));
		}

		writer.close();
	}

}
