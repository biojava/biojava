package org.biojava.nbio.structure;

import java.io.Serializable;

/**
 *  
 * The type of entity (polymer, non-polymer, water, macrolide)
 * as defined in the mmCIF dictionary: <a href="http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v40.dic/Items/_entity.type.html"></a>
 * <p>
 * Entities are of four types:  polymer, non-polymer, macrolide and water.
 * </p>
 * <p>
 * Note that the water category includes only water;  ordered
 * solvent such as sulfate ion or acetone would be described as
 * individual non-polymer entities.
 * It is not clear what macrolides are, but they appear to be supported since mmCIF 4.0.
 * </p>
 *
 * 
 * @author Anthony Bradley
 * @author Jose Duarte
 *
 */
public enum EntityType implements Serializable {

	/**
	 * Polymeric entities: poly-peptides and nucleotide chains
	 */
	POLYMER("polymer"),
	
	/**
	 * Non-polymeric entities: ligands, metal ions, buffer molecules, etc
	 */
	NONPOLYMER("non-polymer"), 
	
	/**
	 * Water
	 */
	WATER("water"),
	
	/**
	 * Macrolide. Supported in mmCIF 4.0 dictionary. Not clear what it refers to.
	 */
	MACROLIDE("macrolide");
	
	private String entityType;

	/**
	 * @param entType the type of the Entity
	 */
	private EntityType(String entType) {
		
		this.setEntityType(entType);
		
	}

	/** 
	 * Returns the type of the Entity as a String
	 *
	 * @return String representation of the entity type.
     */
	public String getEntityType() {
		return entityType;
	}


	private void setEntityType(String entityType) {
		this.entityType = entityType;
	}

	/** 
	 * Creates a new EntityType from a String value.
	 * Returns null if entityType is null or not one of the supported
	 * standard types.
	 *
	 * @param entityType String value , should be one of "polymer","non-polymer","water","macolide"
	 * @return an EntityType object
     */
	public static EntityType entityTypeFromString(String entityType)
	{

		if ( entityType == null)
			return null;

		for(EntityType et : EntityType.values())
		{
			if(entityType.equals(et.entityType))
			{
				return et;
			}
		}
		return null;
	}
	
}
