package org.panda.tfactivity;

/**
 * Gets the change in the expression of a gene.
 * 1: upregulated
 * -1: downregulated
 * 0: no significant change
 *
 * @author Ozgun Babur
 */
public interface DiscreteExpressionProvider
{
	public Integer getChange(String gene);
}
