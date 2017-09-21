package org.panda.tfactivity;

import java.util.*;

/**
 * @author Ozgun Babur
 */
public class RandomizedExpressionProvider implements DiscreteExpressionProvider
{
	DiscreteExpressionProvider provider;
	List<String> genes;
	List<String> shuff;
	Map<String, String> mapping;

	public RandomizedExpressionProvider(DiscreteExpressionProvider provider, Set<String> genesToConsider)
	{
		this.provider = provider;

		genes = new ArrayList<>(genesToConsider);
		shuff = new ArrayList<>(genes);

		mapping = new HashMap<>();
	}

	public void shuffle()
	{
		Collections.shuffle(shuff);

		for (int i = 0; i < genes.size(); i++)
		{
			mapping.put(genes.get(i), shuff.get(i));
		}
	}

	@Override
	public Integer getChange(String gene)
	{
		return provider.getChange(mapping.get(gene));
	}
}
