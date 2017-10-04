package org.panda.tfactivity;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Map;
import java.util.Set;

/**
 * Main class for transcription factor activity prediction.
 *
 * @author Ozgun Babur
 */
public class Main
{
	/**
	 * The name of the parameters file.
	 */
	public static final String PARAMETERS_FILENAME = "parameters.txt";

	/**
	 * The default name of the result file.
	 */
	public static final String DEFAULT_OUTPUT_FILENAME = "TF-activity-results.txt";
	public static final String DEFAULT_GRAPH_FILENAME =  "TF-similarity-graph";

	/**
	 * The directory that contains the parameters file.
	 */
	String directory;

	/**
	 * Name of the output file.
	 */
	String outputFile;

	/**
	 * Name of the output graph (sif) file without .sif extension.
	 */
	String outputGraphFile;

	/**
	 * Type of the expression file.
	 */
	ExpressionResourceType expResType;

	/**
	 * Provider object for the expression in discretized form.
	 */
	DiscreteExpressionProvider expProvider;

	/**
	 * Whether to consider the sign of the network relations.
	 */
	boolean signedAnalysis;

	/**
	 * Iterations for P-value estimation.
	 */
	int iterations;

	/**
	 * The FDR cutoff to use for reporting results.
	 */
	double fdrThr;

	/**
	 * The minimum number of targets of a TF that has an expression data, to consider that TF in the analysis.
	 */
	int minimumTargets;

	/**
	 * Threshold for Jaccard similarity to present on the similarity graph.
	 */
	double jaccardSimilarityThrehsold = 0.2;

	/**
	 * Constructor that sets the working directory and initializes data structures.
	 *
	 * @param directory The directory that contains parameters file and where the output will be generated.
	 */
	public Main(String directory)
	{
		this.directory = directory;
	}

	/**
	 * The method that does the job.
	 *
	 * @throws IOException
	 */
	public void runAnalysis() throws IOException
	{
		System.out.println("directory = " + directory);

		// read the parameters file and set the variables
		readParameters(directory);

		if (outputFile == null) outputFile = directory + "/" + DEFAULT_OUTPUT_FILENAME;
		if (outputGraphFile == null) outputGraphFile= directory + "/" + DEFAULT_GRAPH_FILENAME;

		NetworkLoader nl = new NetworkLoader();

		if (signedAnalysis)
		{
			// load the signed SIF network
			Map<String, Map<String, Integer>> network = nl.loadSignedTCGAConsensus();
			nl.cleanSigned(network, expProvider, minimumTargets);

			// run analysis
			SignedAnalysis sa = new SignedAnalysis(network, expProvider);
			sa.run(outputFile, iterations, fdrThr);
		}
		else
		{
			// load the unsigned network
			Map<String, Set<String>> network = nl.loadUnsigned();
			nl.cleanUnsigned(network, expProvider, minimumTargets);

			// run analysis
			UnsignedAnalysis ua = new UnsignedAnalysis(network, expProvider);
			ua.run(outputFile, fdrThr);
		}

		ResultListToGraph rltg = new ResultListToGraph(outputFile, jaccardSimilarityThrehsold);
		rltg.draw(outputGraphFile);
	}

	/**
	 * reads the parameters file and configures parameters.
	 * @param dir the directory that contains parameters file
	 * @throws IOException
	 */
	void readParameters(String dir) throws IOException
	{
		Files.lines(Paths.get(dir + File.separator + PARAMETERS_FILENAME)).
			filter(l -> !l.startsWith("#")).map(l -> l.split("=")).
			forEach(t ->
			{
				// the token before "=" has to be one of the values in the Parameters enum
				Parameter param = Parameter.findEnum(t[0].trim());

				if (param != null)
				{
					try
					{
						// the specific Parameter enum knows how to configure the Main class using the parameter value
						param.reader.read(t[1].trim(), this);
					}
					catch (IOException e)
					{
						throw new RuntimeException(e);
					}
				}
				else
				{
					System.err.println("Unknown parameter = " + t[0].trim());
				}
			});
	}

	enum Parameter
	{
		OUTPUT_FILE((value, main) -> main.outputFile = value),
		EXPRESSION_RESOURCE_TYPE((value, main) -> main.expResType = ExpressionResourceType.get(value)),
		CONSIDER_EDGE_SIGNS((value, main) -> main.signedAnalysis = Boolean.valueOf(value)),
		RANDOM_ITERATIONS((value, main) -> main.iterations = Integer.valueOf(value)),
		FDR_THRESHOLD((value, main) -> main.fdrThr = Double.valueOf(value)),
		MINIMUM_TARGETS((value, main) -> main.minimumTargets = Integer.valueOf(value)),
		JACCARD_SIMILARITY_THRESHOLD((value, main) -> main.jaccardSimilarityThrehsold= Double.valueOf(value)),

		EXPRESSION_FILE((value, main) ->
		{
			value = main.adjustLocation(value, main.directory);

			switch (main.expResType)
			{
				case TCGA:
				{
					main.expProvider = new TCGAExpressionProvider(value);
					break;
				}
				case Custom:
				{
					main.expProvider = new CustomExpressionProvider(value);
					break;
				}
			}
		}),

		GROUPS_FILE((value, main) ->
		{
			value = main.adjustLocation(value, main.directory);

			switch (main.expResType)
			{
				case TCGA:
				case Custom:
				{
					TwoGroupsSampleNames two = new TwoGroupsSampleNames(value);
					two.filterOutMissingSamples(((TCGAExpressionProvider) main.expProvider).er.getSamples());
					((TCGAExpressionProvider) main.expProvider).setTwo(two);
					break;
				}
				default: throw new RuntimeException("Not implemented yet");
			}
		}),

		EXPRESSION_PVAL_THRESHOLD((value, main) ->
		{
			switch (main.expResType)
			{
				case TCGA:
				{
					((TCGAExpressionProvider) main.expProvider).setPvalThr(Double.valueOf(value));
					break;
				}
				case Custom: throw new RuntimeException("Not implemented yet");
			}
		}),

		EXPRESSION_FDR_THRESHOLD((value, main) ->
		{
			switch (main.expResType)
			{
				case TCGA:
				case Custom:
				{
					((TCGAExpressionProvider) main.expProvider).setFDRThr(Double.valueOf(value));
					break;
				}
				default: throw new RuntimeException("Not implemented yet");
			}
		}),
		;

		ParameterReader reader;

		Parameter(ParameterReader reader)
		{
			this.reader = reader;
		}

		String getText()
		{
			return toString().toLowerCase().replaceAll("_", "-");
		}

		static Parameter findEnum(String text)
		{
			for (Parameter parameter : values())
			{
				if (parameter.getText().equals(text)) return parameter;
			}
			return null;
		}
	}

	private String adjustLocation(String file, String baseDir)
	{
		if (!file.startsWith("/")) return baseDir + "/" + file;
		return file;
	}

	interface ParameterReader
	{
		void read(String value, Main main) throws IOException;
	}

	public static void main(String[] args) throws IOException
	{
		Main main = new Main(args[0]);
		main.runAnalysis();
	}
}
