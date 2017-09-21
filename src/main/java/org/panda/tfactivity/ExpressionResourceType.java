package org.panda.tfactivity;

/**
 * @author Ozgun Babur
 */
public enum ExpressionResourceType
{
	TCGA,
	Custom;

	public static ExpressionResourceType get(String val)
	{
		for (ExpressionResourceType type : values())
		{
			if (type.toString().equals(val)) return type;
		}
		return null;
	}
}
