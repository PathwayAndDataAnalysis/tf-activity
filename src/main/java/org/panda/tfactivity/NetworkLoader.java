package org.panda.tfactivity;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.resource.network.*;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.graph.DirectedGraph;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Loads Pathway Commons SIF.
 *
 * @author Ozgun Babur
 */
public class NetworkLoader
{
	/**
	 * Loads the network.
	 * @return ordered sets of graphs
	 */
	public Map<String, Map<String, Integer>> loadSigned()
	{
		Map<String, Map<String, Integer>> map = new HashMap<>();

		populateSigned(map, SignedPC.get().getGraph(SignedType.UPREGULATES_EXPRESSION), 1);
		populateSigned(map, SignedPC.get().getGraph(SignedType.DOWNREGULATES_EXPRESSION), -1);
		populateSigned(map, TRRUST.get().getPositiveGraph(), 1);
		populateSigned(map, TRRUST.get().getNegativeGraph(), -1);
		populateSigned(map, TFactS.get().getPositiveGraph(), 1);
		populateSigned(map, TFactS.get().getNegativeGraph(), -1);

		// remove conflicts

		for (String tf : map.keySet())
		{
			Set<String> remove = map.get(tf).keySet().stream().filter(target -> map.get(tf).get(target) == 0)
				.collect(Collectors.toSet());

			remove.forEach(map.get(tf)::remove);
		}

		return map;
	}

	private void populateSigned(Map<String, Map<String, Integer>> map, DirectedGraph graph, int rel)
	{
		for (String tf : graph.getOneSideSymbols(true))
		{
			if (!map.containsKey(tf)) map.put(tf, new HashMap<>());

			for (String target : graph.getDownstream(tf))
			{
				if (map.get(tf).containsKey(target))
				{
					// if there is a conflict, mark it with zero
					if (map.get(tf).get(target) * rel == -1)
					{
						map.get(tf).put(target, 0);
					}
				}
				else
				{
					map.get(tf).put(target, rel);
				}
			}
		}
	}

	public Map<String, Set<String>> loadUnsigned()
	{
		Map<String, Set<String>> map = new HashMap<>();
		populateUnsigned(map, (DirectedGraph) PathwayCommons.get().getGraph(SIFEnum.CONTROLS_EXPRESSION_OF));
		populateUnsigned(map, TRRUST.get().getUnsignedGraph());
		populateUnsigned(map, TFactS.get().getUnsignedGraph());
		return map;
	}

	private void populateUnsigned(Map<String,Set<String>> map, DirectedGraph graph)
	{
		for (String tf : graph.getOneSideSymbols(true))
		{
			if (!map.containsKey(tf)) map.put(tf, new HashSet<>());
			map.get(tf).addAll(graph.getDownstream(tf));
		}
	}

	public void cleanSigned(Map<String, Map<String, Integer>> network, DiscreteExpressionProvider prov, int minTargets)
	{
		// remove targets with no expression data
		for (String tf : network.keySet())
		{
			new HashSet<>(network.get(tf).keySet()).stream().filter(target -> prov.getChange(target) == null)
				.forEach(target -> network.get(tf).remove(target));
		}

		// remove TFs with low number of targets
		new HashSet<>(network.keySet()).stream()
			.filter(tf -> network.get(tf).size() < minTargets)
			.forEach(network::remove);
	}

	public void cleanUnsigned(Map<String, Set<String>> network, DiscreteExpressionProvider prov, int minTargets)
	{
		// remove targets with no expression data
		for (String tf : network.keySet())
		{
			new HashSet<>(network.get(tf)).stream().filter(target -> prov.getChange(target) == null)
				.forEach(target -> network.get(tf).remove(target));
		}

		// remove TFs with low number of targets
		new HashSet<>(network.keySet()).stream()
			.filter(tf -> network.get(tf).size() < minTargets)
			.forEach(network::remove);
	}
}
