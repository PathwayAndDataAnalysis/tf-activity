package org.panda.tfactivity;

import org.panda.utility.Progress;
import org.panda.utility.statistics.FDR;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class SignedAnalysis
{
	/**
	 * Signed transcriptional network.
	 */
	Map<String, Map<String, Integer>> network;

	/**
	 * Expression data provider.
	 */
	DiscreteExpressionProvider expProv;

	public SignedAnalysis(Map<String, Map<String, Integer>> network, DiscreteExpressionProvider expProv)
	{
		this.network = network;
		this.expProv = expProv;
	}

	public void run(String outFile, int iterations, double fdrThr) throws IOException
	{
		Map<String, int[]> counts = getActivityEvidenceCounts(network, expProv);
		Map<String, double[]> pValues = getActivityPValues(network, expProv, counts, iterations);
		writeResults(outFile, network, expProv, pValues, fdrThr);
	}

	private Map<String, double[]> getActivityPValues(Map<String, Map<String, Integer>> network,
		DiscreteExpressionProvider expProv, Map<String, int[]> actualCounts, int iterations)
	{
		Set<String> genesToShuffle = network.values().stream().map(Map::keySet)
			.flatMap(Collection::stream).collect(Collectors.toSet());

		RandomizedExpressionProvider randProv = new RandomizedExpressionProvider(expProv, genesToShuffle);

		Map<String, int[]> betterCnts = new HashMap<>();
		for (String tf : actualCounts.keySet())
		{
			betterCnts.put(tf, new int[]{0, 0});
		}

		Progress p = new Progress(iterations, "Calculating p-values");
		for (int i = 0; i < iterations; i++)
		{
			randProv.shuffle();

			Map<String, int[]> cnts = getActivityEvidenceCounts(network, randProv);

			for (String tf : cnts.keySet())
			{
				int[] rand = cnts.get(tf);
				int[] act = actualCounts.get(tf);
				int[] bett = betterCnts.get(tf);

				if (rand[0] >= act[0]) bett[0]++;
				if (rand[1] >= act[1]) bett[1]++;
			}
			p.tick();
		}

		Map<String, double[]> pvalMap = new HashMap<>();
		for (String tf : betterCnts.keySet())
		{
			pvalMap.put(tf, new double[]{
				betterCnts.get(tf)[0] / (double) iterations, betterCnts.get(tf)[1] / (double) iterations});
		}

		return pvalMap;
	}

	private Map<String, int[]> getActivityEvidenceCounts(Map<String, Map<String, Integer>> network, 
		DiscreteExpressionProvider expProv)
	{
		Map<String, int[]> evidenceMap = new HashMap<>();

		for (String tf : network.keySet())
		{
			int[] cnt = new int[]{0, 0};

			Map<String, Integer> map = network.get(tf);

			for (String target : map.keySet())
			{
				int change = expProv.getChange(target);

				if (change != 0)
				{
					int edgeSign = map.get(target);
					cnt[edgeSign * change > 0 ? 0 : 1]++;
				}
			}

			evidenceMap.put(tf, cnt);
		}

		return evidenceMap;
	}

	private Map<String, List<String>[]> getActivitySupporterGenes(Map<String, Map<String, Integer>> network,
		DiscreteExpressionProvider expProv)
	{
		Map<String, List<String>[]> evidenceMap = new HashMap<>();

		for (String tf : network.keySet())
		{
			List<String>[] lists = new List[]{new ArrayList<String>(), new ArrayList<String>()};

			Map<String, Integer> map = network.get(tf);

			for (String target : map.keySet())
			{
				int change = expProv.getChange(target);

				if (change != 0)
				{
					int edgeSign = map.get(target);
					lists[edgeSign * change > 0 ? 0 : 1].add(target);
				}
			}

			Collections.sort(lists[0]);
			Collections.sort(lists[1]);

			evidenceMap.put(tf, lists);
		}

		return evidenceMap;
	}

	private void writeResults(String file, Map<String, Map<String, Integer>> network,
		DiscreteExpressionProvider expProv, Map<String, double[]> pvals, double fdrThr) throws IOException
	{
		Map<String, Double> pMap = new HashMap<>();
		for (String tf : pvals.keySet())
		{
			pMap.put(tf + "-act", pvals.get(tf)[0]);
			pMap.put(tf + "-inh", pvals.get(tf)[1]);
		}

		Map<String, Double> qMap = FDR.getQVals(pMap, null);

		List<String> select = FDR.select(pMap, null, fdrThr);

		Map<String, List<String>[]> support = getActivitySupporterGenes(network, expProv);

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(file));
		writer.write("TF\tChange\tP-val\tQ-val\tSupporting targets");

		for (String s : select)
		{
			String tf = s.substring(0, s.lastIndexOf("-"));
			String change = s.substring(s.lastIndexOf("-") + 1).equals("act") ? "activated" : "inhibited";

			writer.write("\n" + tf + "\t" + change + "\t" + pMap.get(s) + "\t" + qMap.get(s) + "\t" +
				support.get(tf)[change.startsWith("a") ? 0 : 1]);
		}

		writer.close();
	}

}
