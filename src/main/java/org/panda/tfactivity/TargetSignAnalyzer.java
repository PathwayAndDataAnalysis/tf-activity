package org.panda.tfactivity;

import org.panda.resource.tcga.ExpressionReader;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.Kronometre;
import org.panda.utility.Tuple;
import org.panda.utility.statistics.Binomial;
import org.panda.utility.statistics.Correlation;
import org.panda.utility.statistics.FDR;
import org.panda.utility.statistics.Summary;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Infers TF --> target edge signs from TCGA expressions.
 *
 * @author Ozgun Babur
 */
public class TargetSignAnalyzer
{
	public static void main(String[] args) throws IOException
	{
		Kronometre k = new Kronometre();
		String dir = "/home/babur/Documents/PC/";
		TargetSignAnalyzer tsa = new TargetSignAnalyzer();
		tsa.run(0.05, 0, dir + "SignedByTCGAConsensus.sif", dir + "SignedByTCGAConsensusFiltered.sif", dir +
			"PC-TCGA-expression-agreement-analysis.txt", dir);
		k.print();
	}
	void run(double fdrThr, double expStdevThr, String unfilteredOutFile, String filteredOutFile,
		String agreementAnalysisFile, String dir) throws IOException
	{
		NetworkLoader nl = new NetworkLoader();
		Map<String, Set<String>> unsignedMap = nl.loadUnsigned();
		Map<String, Map<String, Integer>> signedMap = nl.loadSigned();
//		ensureUnsignedContainsAllSigned(unsignedMap, signedMap);

		Set<String> genes = unsignedMap.values().stream().flatMap(Collection::stream).collect(Collectors.toSet());
		genes.addAll(unsignedMap.keySet());

		Map<String, Map<String, double[]>> expsMap = getExpressions(genes, expStdevThr);
		System.out.println("Loaded " + expsMap.size() + " expression sets.");

		if (Math.random() < 0)
		{
			TFExpToActMapper mapper = new TFExpToActMapper(signedMap, expsMap, dir + "/MonotCheck");
			mapper.run();

			return;
		}

		// Count positive and negative correlations
		Map<String, Map<String, Pair>> cntMap = countCorrelations(unsignedMap, expsMap, fdrThr, expStdevThr);
		System.out.println("Counted all correlations");

		// Decide consensus regulation.
		assignDecisions(cntMap, fdrThr);
		System.out.println("Decided consensus.");

		// Write results
		writeSignificantResults(cntMap, unfilteredOutFile);
		System.out.println("Wrote unfiltered results.");

		Set<String> behavedFactors = checkAgreementWithSignedPC(cntMap, signedMap, fdrThr, agreementAnalysisFile);
		System.out.println("Checked agreement with SignedPC.");

		Set<String> remove = cntMap.keySet().stream().filter(f -> !behavedFactors.contains(f))
			.collect(Collectors.toSet());

		remove.forEach(cntMap::remove);

		writeSignificantResults(cntMap, filteredOutFile);
		System.out.println("Wrote filtered.");
	}

	/**
	 * Returns the names of factors whose majority of biased targets are biased towards the expected direction.
	 */
	private Set<String> checkAgreementWithSignedPC(Map<String, Map<String, Pair>> cntMap,
		Map<String, Map<String, Integer>> signedMap, double fdrThr, String outFile) throws IOException
	{
		Map<String, Integer> fAgr = new HashMap<>();
		Map<String, Integer> fDis = new HashMap<>();
		int totAgr = 0;
		int totDis = 0;

		for (String factor : cntMap.keySet())
		{
			if (signedMap.containsKey(factor))
			{
				if (!fAgr.containsKey(factor)) fAgr.put(factor, 0);
				if (!fDis.containsKey(factor)) fDis.put(factor, 0);

				Map<String, Integer> fSPC = signedMap.get(factor);

				for (Pair pair : cntMap.get(factor).values())
				{
					if (pair.consensus != null && fSPC.containsKey(pair.target))
					{
						if (pair.consensus.equals(fSPC.get(pair.target)))
						{
							fAgr.put(factor, fAgr.get(factor) + 1);
							totAgr++;
						}
						else
						{
							fDis.put(factor, fDis.get(factor) + 1);
							totDis++;
						}
					}
				}
			}
		}
		
		// Detect significantly deviated
		
		Map<String, Double> biasPval = new HashMap<>();
		Map<String, Double> biasLimit = new HashMap<>();
		Set<String> posBias = new HashSet<>();
		for (String factor : fAgr.keySet())
		{
			int agree = fAgr.get(factor);
			int confl = fDis.get(factor);
			int total = agree + confl;
			if (total > 0)
			{
				if (agree > confl)
				{
					posBias.add(factor);
				}

				biasPval.put(factor, Binomial.getPval(agree, confl));
				biasLimit.put(factor, Binomial.getPval(total, 0));
			}
		}
		Set<String> select = new HashSet<>(FDR.select(biasPval, biasLimit, fdrThr));

		// Write down agree/conflict numbers
		
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));

		System.out.println("totAgg = " + totAgr);
		System.out.println("totDis = " + totDis);
		writer.write("Total common targets = " + (totAgr + totDis));
		writer.write("\nTotal agreement = " + (totAgr / (double) (totAgr + totDis)));

		writer.write("\n\nFactor\tCommon targets\tAgreement\tSignificant");
		for (String factor : fAgr.keySet())
		{
			int total = fAgr.get(factor) + fDis.get(factor);
			if (total > 0)
			{
				writer.write("\n" + factor + "\t" + total + "\t" + (fAgr.get(factor) / (double) total));
				if (select.contains(factor)) writer.write("\tX");
			}
		}
		writer.close();

		Set<String> negSignif = CollectionUtil.diff(select, posBias);
		System.out.println("Genes that are significantly reverse correlated = " + negSignif);

		select.retainAll(posBias);
		return select;
	}

	private void writeSignificantResults(Map<String, Map<String, Pair>> cntMap, String filename) throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(filename));

		for (String factor : cntMap.keySet())
		{
			Map<String, Pair> fMap = cntMap.get(factor);

			for (String target : fMap.keySet())
			{
				Pair pair = fMap.get(target);

				if (pair.consensus != null)
				{
					writer.write(factor + "\t" + (pair.consensus > 0 ? "up" : "down") + "regulates-expression\t" +
						target + "\n");
				}
			}
		}

		writer.close();
	}

	void assignDecisions(Map<String, Map<String, Pair>> cntMap, double fdrThr)
	{
		for (String factor : cntMap.keySet())
		{
			Map<String, Pair> fCnt = cntMap.get(factor);

			Map<String, Double> pvals = new HashMap<>();
			Map<String, Double> limits = new HashMap<>();

			for (String target : fCnt.keySet())
			{
				Pair pair = fCnt.get(target);
				pvals.put(target, pair.getBinomialPvalue());
				limits.put(target, pair.getMaxPossibleBinomPval());
			}

			double thr = FDR.getPValueThreshold(pvals, limits, fdrThr);
			pvals.keySet().stream().filter(t -> pvals.get(t) < thr).forEach(t ->
			{
				Pair pair = fCnt.get(t);
				pair.consensus = (pair.pos > pair.neg) ? 1 : -1;
			});
		}
	}

	private Map<String, Map<String, Pair>> countCorrelations(Map<String, Set<String>> unsignedMap,
		Map<String, Map<String, double[]>> expsMap, double fdrThr, double expStdevThr)
	{
		Map<String, Map<String, Pair>> cntMap = new HashMap<>();

		for (String factor : unsignedMap.keySet())
		{
			for (Map<String, double[]> exps : expsMap.values())
			{
				Map<String, Integer> dExp = getDiscretizedExpressions(factor, unsignedMap.get(factor), exps, fdrThr,
					expStdevThr);

				if (dExp != null && !dExp.isEmpty())
				{
					if (!cntMap.containsKey(factor)) cntMap.put(factor, new HashMap<>());

					Map<String, Pair> fMap = cntMap.get(factor);

					for (String target : dExp.keySet())
					{
						Pair pair = fMap.get(target);
						if (pair == null)
						{
							pair = new Pair(factor, target);
							fMap.put(target, pair);
						}

						if (dExp.get(target) == 1) pair.pos++;
						else
						{
							assert dExp.get(target) == -1;
							pair.neg++;
						}
					}
				}
			}
		}
		return cntMap;
	}

	Map<String, Integer> getDiscretizedExpressions(String factor, Set<String> targets, Map<String, double[]> exps,
		double fdrThr, double expStdevThr)
	{
		Map<String, Double> pvals = new HashMap<>();
		Map<String, Integer> directions = new HashMap<>();

		double[] fval = exps.get(factor);

		if (fval != null)
		{
			if (Summary.stdev(fval) < expStdevThr) return null;

			for (String target : targets)
			{
				double[] tval = exps.get(target);
				if (tval != null)
				{
					double[][] v = ArrayUtil.trimNaNs(fval, tval);
					Tuple cor = Correlation.pearson(v[0], v[1]);
					pvals.put(target, cor.p);
					if (cor.v != 0) directions.put(target, (int) Math.signum(cor.v));
				}
			}

			if (!pvals.isEmpty())
			{
				double thr = FDR.getPValueThreshold(pvals, null, fdrThr);
				pvals.keySet().stream().filter(t -> pvals.get(t) > thr).forEach(directions::remove);
			}
		}
		return directions;
	}

	Map<String, Map<String, double[]>> getExpressions(Set<String> genes, double stdevThr)
	{
		Map<String, Map<String, double[]>> map = new HashMap<>();

		File dir = new File("/home/babur/Documents/TCGA");
		Arrays.stream(dir.listFiles()).filter(f -> !containsUnwanted(f.getName())).forEach(f -> {
			try
			{
				ExpressionReader er = new ExpressionReader(f.getPath() + "/expression.txt");
				String[] samples = er.getSamples().toArray(new String[0]);
				Map<String, double[]> exps = new HashMap<>();
				for (String gene : genes)
				{
					double[] exp = er.getGeneAlterationArray(gene, samples);
					if (exp != null)
					{
//						double stdev = Summary.stdev(exp);
//						if (stdev >= stdevThr)
							exps.put(gene, exp);
					}
				}
				map.put(f.getName(), exps);
			}
			catch (FileNotFoundException e)
			{
				throw new RuntimeException(e);
			}
		});

		return map;
	}

	boolean containsUnwanted(String file)
	{
		return file.contains("GBMLGG") || file.contains("COADREAD") || file.contains("KIPAN") ||
			file.contains("PanCan") || file.contains("KIPAN") || file.contains("STAD");
	}

	void ensureUnsignedContainsAllSigned(Map<String, Set<String>> unsignedMap,
		Map<String, Map<String, Integer>> signedMap)
	{
		for (String tf : signedMap.keySet())
		{
			if (!unsignedMap.containsKey(tf)) throw new RuntimeException("Unsigned is missing tf: " + tf);

			if (!unsignedMap.get(tf).containsAll(signedMap.get(tf).keySet()))
			{
				Set<String> diff = CollectionUtil.diff(signedMap.get(tf).keySet(), unsignedMap.get(tf));
				throw new RuntimeException("Missing " + tf + " targets in unsigned = " + diff);
			}
		}
	}

	class Pair
	{
		String factor;
		String target;
		int pos;
		int neg;
		Integer consensus;

		public Pair(String factor, String target)
		{
			this.factor = factor;
			this.target = target;
			pos = 0;
			neg = 0;
		}

		double getBinomialPvalue()
		{
			return Binomial.getPval(pos, neg);
		}

		double getMaxPossibleBinomPval()
		{
			return Binomial.getPval(pos + neg, 0);
		}
	}
}
