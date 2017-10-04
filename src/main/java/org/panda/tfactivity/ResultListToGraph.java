package org.panda.tfactivity;

import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.ValToColor;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class ResultListToGraph
{
	String inFile;
	double minimumJaccardSimilarity;

	double mostSignifPval = 1E-5;
	int edgeWeights = 2;
	Color maxUp = new Color(255, 100, 100);
	Color maxDw = new Color(100, 100, 255);
	Color maxSig = new Color(50, 200, 50);
	Color neutr = Color.WHITE;

	int tfInd;
	int pInd;
	int actInd;
	int targetInd;

	Map<String, Row> data;

	public ResultListToGraph(String inFile, double minimumJaccardSimilarity) throws IOException
	{
		this.inFile = inFile;
		this.minimumJaccardSimilarity = minimumJaccardSimilarity;
		load();
	}

	void load() throws IOException
	{
		String[] header = Files.lines(Paths.get(inFile)).findFirst().get().split("\t");
		tfInd = ArrayUtil.indexOf(header, "TF");
		pInd = ArrayUtil.indexOf(header, "P-val");
		actInd = ArrayUtil.indexOf(header, "Change");
		targetInd = ArrayUtil.indexOf(header, "Supporting targets");

		data = new HashMap<>();
		Scanner sc = new Scanner(new File(inFile));
		sc.nextLine();
		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			Row row = new Row(line);

			if (data.containsKey(row.id))
			{
				Row other = data.remove(row.id);
				other.attachActivityToID();
				data.put(other.id, other);
				row.attachActivityToID();
			}
			data.put(row.id, row);
		}
	}

	public void draw(String outFileWoExt) throws IOException
	{
		ValToColor edgeCol = new ValToColor(new double[]{0, 1}, new Color[]{Color.WHITE, Color.BLACK});
		double maxScore = -Math.log(mostSignifPval);
		ValToColor nodeCol = data.values().iterator().next().activity != null ?
			new ValToColor(new double[]{-maxScore, 0, maxScore}, new Color[]{maxDw, neutr, maxUp}) :
			new ValToColor(new double[]{0, maxScore}, new Color[]{neutr, maxSig});
		String edgeType = "correlates-with";

		BufferedWriter sifWriter = Files.newBufferedWriter(Paths.get(outFileWoExt + ".sif"));
		BufferedWriter fmtWriter = Files.newBufferedWriter(Paths.get(outFileWoExt + ".format"));
		fmtWriter.write("edge\tall-edges\twidth\t" + edgeWeights);

		for (String id1 : data.keySet())
		{
			Row tf1 = data.get(id1);
			for (String id2 : data.keySet())
			{
				if (id1.compareTo(id2) < 0)
				{
					Row tf2 = data.get(id2);

					double sim = tf1.getSimilarity(tf2);

					if (sim >= minimumJaccardSimilarity)
					{
						sifWriter.write(id1 + "\t" + edgeType + "\t" + id2 + "\n");
						fmtWriter.write("\nedge\t" + id1 + " " + edgeType + " " + id2 + "\tcolor\t" +
							edgeCol.getColorInString(sim));
					}
				}
			}
		}
		data.keySet().forEach(id ->
		{
			FileUtil.writeln(id, sifWriter);
			FileUtil.lnwrite("node\t" + id + "\tcolor\t" + nodeCol.getColorInString(
				-Math.log(data.get(id).p) * (data.get(id).activity != null ? data.get(id).activity : 1)),
				fmtWriter);
		});

		sifWriter.close();
		fmtWriter.close();
	}

	class Row
	{
		String id;
		String tf;
		double p;
		Integer activity;
		Set<String> targets;

		public Row(String s)
		{
			String[] t = s.split("\t");
			tf = t[tfInd];
			id = tf;
			p = Double.valueOf(t[pInd]);
			if (actInd >= 0) activity = t[actInd].startsWith("a") ? 1 : -1;
			t[targetInd] = t[targetInd].substring(1, t[targetInd].length() - 1);
			targets = new HashSet<>(Arrays.asList(t[targetInd].split(", ")));
		}

		void attachActivityToID()
		{
			id += "-" + (activity == 1 ? "a" : "i");
		}

		double getSimilarity(Row other)
		{
			return CollectionUtil.getJaccardSimilarity(targets, other.targets);
		}
	}

	public static void main(String[] args) throws IOException
	{
		String dir = "/home/babur/Documents/Analyses/TF-activity/TCGA-OV/Immunoreactive/Signed/";
		ResultListToGraph rltg = new ResultListToGraph(dir + "TF-activity-results.txt", 0.2);
		rltg.draw(dir + "result-graph");
	}
}
