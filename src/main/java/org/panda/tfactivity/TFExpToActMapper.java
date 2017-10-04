package org.panda.tfactivity;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class TFExpToActMapper
{
	Map<String, Map<String, Integer>> signedMap;
	Map<String, Map<String, double[]>> expsMap;
	Map<String, Map<String, double[]>> ranksMap;
	Map<String, Map<String, double[]>> ranksRevMap;

	String outDir;

	public TFExpToActMapper(Map<String, Map<String, Integer>> signedMap, Map<String, Map<String, double[]>> expsMap,
		String outDir)
	{
		this.signedMap = signedMap;
		this.expsMap = expsMap;
		this.outDir = outDir;
		ranksMap = new HashMap<>();
		ranksRevMap = new HashMap<>();
	}

	public void run() throws IOException
	{
		for (String factor : signedMap.keySet())
		{
			Map<String, Integer> fMap = signedMap.get(factor);

			if (fMap.size() > 5)
			{
				for (String study : expsMap.keySet())
				{
					double[] fR = getRanks(study, factor, true);

					if (fR == null || fR.length < 10) continue;

					double[] ac = new double[fR.length];

					for (int i = 0; i < fR.length; i++)
					{
						double sum = 0;
						double cnt = 0;

						for (String target : fMap.keySet())
						{
							boolean sign = fMap.get(target) == 1;
							double[] tR = getRanks(study, target, sign);
							if (tR != null && !Double.isNaN(tR[i]))
							{
								sum += tR[i];
								cnt++;
							}
						}

						ac[i] = sum / cnt;
					}

					double[][] arr = applySlidingWindow(fR, ac, 0.1);
					write(study, factor, arr[0], arr[1]);
				}
			}
		}
	}

	void write(String study, String factor, double[] fRank, double[] fAc) throws IOException
	{
		File dir = new File(outDir + "/" + factor);
		dir.mkdirs();

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(dir.getPath() + "/" + study + ".txt"));
		writer.write(study + "\t" + factor + "\nFactor rank\tTarget expression");

		for (int i = 0; i < fRank.length; i++)
		{
			writer.write("\n" + fRank[i] + "\t" + fAc[i]);
		}

		writer.close();
	}

	double[] getRanks(String study, String gene, boolean forward)
	{
		if (forward && !ranksMap.containsKey(study)) ranksMap.put(study, new HashMap<>());
		if (!forward && !ranksRevMap.containsKey(study)) ranksRevMap.put(study, new HashMap<>());

		if (forward && !ranksMap.get(study).containsKey(gene) ||
			!forward && !ranksRevMap.get(study).containsKey(gene))
		{
			double[] exp = expsMap.get(study).get(gene);

			if (exp != null)
			{

				Map<String, Map<String, double[]>> map = forward ? ranksMap : ranksRevMap;
				map.get(study).put(gene, getRanks(exp, forward));
			}
		}

		return forward ? ranksMap.get(study).get(gene) : ranksRevMap.get(study).get(gene);
	}

	double[] getRanks(double[] vals, boolean forward)
	{
		List<Node> nodes = new ArrayList<>(vals.length);
		for (int i = 0; i < vals.length; i++)
		{
			nodes.add(new Node(vals[i], i, forward));
		}

		Collections.sort(nodes);

		for (int i = 0; i < nodes.size(); i++)
		{
			Node node = nodes.get(i);

			if (node.val.isNaN())
			{
				node.rank = Double.NaN;
			}
			else if (i == 0)
			{
				node.rank = i;
			}
			else
			{
				Node prev = nodes.get(i - 1);

				if (prev.val.equals(node.val)) node.rank = prev.rank;
				else node.rank = i;
			}
		}

		double[] r = new double[vals.length];

		for (Node node : nodes)
		{
			r[node.index] = node.rank;
		}

		return r;
	}

	double[][] applySlidingWindow(double[] fR, double[] fA, double winRatio)
	{
		// sort arrays using fA

		List<Tuple> list = new ArrayList<>();
		for (int i = 0; i < fR.length; i++)
		{
			list.add(new Tuple(fR[i], fA[i]));
		}
		Collections.sort(list);
		fR = new double[fR.length];
		fA = new double[fA.length];
		for (int i = 0; i < fR.length; i++)
		{
			fR[i] = list.get(i).v1;
			fA[i] = list.get(i).v2;
		}

		int window = (int) (fR.length * winRatio);

		double[] fRW = new double[fR.length - window + 1];
		double[] fAW = new double[fRW.length];

		for (int i = 0; i < fRW.length; i++)
		{
			fRW[i] = fR[i + window - 1];

			for (int j = i; j < i + window; j++)
			{
				fAW[i] += fA[j];
			}
		}

		for (int i = 0; i < fRW.length; i++)
		{
			fAW[i] /= window;
		}

		return new double[][]{fRW, fAW};
	}

	class Tuple implements Comparable
	{
		Double v1;
		Double v2;

		public Tuple(Double v1, Double v2)
		{
			this.v1 = v1;
			this.v2 = v2;
		}

		@Override
		public int compareTo(Object o)
		{
			return v1.compareTo(((Tuple) o).v1);
		}
	}

	class Node implements Comparable
	{
		Double val;
		int index;
		double rank;
		boolean forward;

		public Node(double val, int index, boolean forward)
		{
			this.val = val;
			this.index = index;
			this.forward = forward;
		}

		@Override
		public int compareTo(Object o)
		{
			if (forward) return val.compareTo(((Node) o).val);

			if (val.isNaN() && !((Node) o).val.isNaN()) return -1;
			if (!val.isNaN() && ((Node) o).val.isNaN()) return 1;
			return ((Node) o).val.compareTo(val);
		}
	}
}
