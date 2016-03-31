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

	/**
	 * @param entType the type of the Entity
	 */
	private EntityType(String entType) {
		
		this.setEntityType(entType);
		
	}

	/** Returns the type of the Entity as a String
	 *
	 * @return String representation of the entity type.
     */
	public String getEntityType() {
		return entityType;
	}


	private void setEntityType(String entityType) {
		this.entityType = entityType;
	}

	/** Creates a new EntityType from a String value.
	 *  Returns null if entityType is null or not one of the supported
	 *  standard types.
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
