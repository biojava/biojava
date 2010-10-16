/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */

package org.biojava.directory;

import java.util.Collections;
import java.util.Map;

import org.biojava.utils.OverlayMap;

/**
 * <p>The BioDirectory Registry is a simple system for specifying
 * where to find services which provide sequence databases. A client
 * should look at the following places in order to find this file:</p>
 *
 * <pre>
 *  $HOME/.bioinformatics/seqdatabase.ini
 *  /etc/bioinformatics/seqdatabase.ini
 *  http://www.open-bio.org/registry/seqdatabase.ini
 * </pre>
 *
 * <p>The file is a simple stanza format</p>
 *
 * <pre>
 *  [database-name]
 *  tag=value
 *  tag=value
 *
 *  [database-name]
 *  tag=value
 *  tag=value
 * </pre>
 * 
 * <p>where each stanza starts with the declaration of the database
 * being in square brackets and following that one line tag=value tag
 * value formats.</p>
 *
 * <p>Database-name stanzas can be repeated, in which case the client
 * should try each service in turn (starting at the first one).</p>
 *
 * <p>The options under each stanza must have two non-optional
 * tag=value lines being</p>
 *
 * <pre>
 *  protocol=&lt;protocol-type&gt;
 *  location=&lt;location-string&gt;
 * </pre>
 *
 * <p>'protocol' currently can be one of</p>
 *
 * <ul>
 *  <li>flat</li>
 *  <li>biofetch</li>
 *  <li>biosql</li>
 * </ul>
 *
 * <p>'location' is a string specific to the protocol. Any number of
 * additional tag values are allowed in the stanza which should be
 * passed to the appropiate constructor of the protocol to
 * interpret. Some protocols might insist on other mandatory tags.</p>
 *
 * @author Brian Gilman
 * @author Keith James
 * @author Matthew Pocock
 * @version $Revision$
 */
public interface RegistryConfiguration {
  /**
   * <code>getConfiguration</code> returns a mapping of registry
   * database names to collections of tag-value pairs.
   *
   * @return a <code>Map</code>.
   *
   * @exception RegistryException if an error occurs.
   */
  public Map getConfiguration() throws RegistryException;
  
  /**
   * <code>getConfigLocator</code> returns a locator for the
   * configuration.
   *
   * @return a <code>String</code>.
   */
  public String getConfigLocator();
  
  /**
   * A simple implementation of RegistryConfiguration backed by a Map.
   *
   * @author Brian Gilman
   * @author Matthew Pocock
   */
  public static class Impl implements RegistryConfiguration {
    private String configFileLocation = null;
    private Map config = null;
    
    public Impl(String configFileLocation, Map config){
      this.configFileLocation = configFileLocation;
      this.config = config;
    }

    public Map getConfiguration() {
      return config;
    }

    public String getConfigLocator() {
      return configFileLocation;
    }
  }
  
  /**
   * A RegistryConfiguration that allows you to treat other
   * configurations as providing important or default configuration
   * information.
   *
   * @author Matthew Pocock
   */
  public static class Composite
  implements RegistryConfiguration {
    private String configLocator;
    private Map config;

    public Composite() {
    }

    public Map getConfiguration() {
      if(config == null) {
        return Collections.EMPTY_MAP;
      } else {
        return config;
      }
    }

    public String getConfigLocator() {
      return configLocator;
    }

    /**
     * Add a configuration as the most authoritative place to look.
     * During future lookups with this context, values in newConfig
     * will take precedence over values in the previously existing
     * configuration.
     *
     * @param newConfig the RegistryConfiguration to add as most
     * important
     */
    public void addTopConfig(RegistryConfiguration newConfig)
    throws RegistryException {
      Map cfg = newConfig.getConfiguration();
      if(config == null) {
        config = cfg;
        configLocator = newConfig.getConfigLocator();
      } else {
        config = new OverlayMap(config, cfg);
        configLocator = newConfig.getConfigLocator() + "::" + configLocator;
      }
    }

    /**
     * Add a configuration as the most default place to look. During
     * future lookups with this context, values in newConfig will be
     * used as default values only if the lookup would return nothing
     * in the previously existing configuration.
     *
     * @param newConfig the RegistryConfiguration to add as the
     * default
     */
    public void addBottomConfig(RegistryConfiguration newConfig)
    throws RegistryException {
      Map cfg = newConfig.getConfiguration();
      if(config == null) {
        config = cfg;
        configLocator = newConfig.getConfigLocator();
      } else {
        config = new OverlayMap(cfg, config);
        configLocator = configLocator + "::" + newConfig.getConfigLocator();
      }
    }
  }
}
