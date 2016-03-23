package org.biojava.nbio.structure;

import java.io.Serializable;

/**
 * This is the type of entity.
 * @author Anthony Bradley
 *
 */
public enum EntityType implements Serializable {

	/**
	 * The type of entity (polymer, non-polymer, water, macrolide)
	 * Defined here:http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v40.dic/Items/_entity.type.html
	 * """Entities are of three types:  polymer, non-polymer and water.
               Note that the water category includes only water;  ordered
               solvent such as sulfate ion or acetone would be described as
               individual non-polymer entities."""
               NOT SURE WHAT MACROLIDES ARE BUT THEY APPEAR TO BE SUPPORTED.....
	 */
	POLYMER("polymer"),
	NONPOLYMER("non-polymer"), 
	WATER("water"),
	MACROLIDE("macrolide");
	
	private String entityType;
	
	private EntityType(String entType) {
		
		this.setEntityType(entType);
		
	}

	public String getEntityType() {
		return entityType;
	}

	public void setEntityType(String entityType) {
		this.entityType = entityType;
	}
	
	public static EntityType entityTypeFromString(String entityType)
	{
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
