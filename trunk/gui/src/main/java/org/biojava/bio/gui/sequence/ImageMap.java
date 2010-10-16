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

package org.biojava.bio.gui.sequence;

import java.io.Serializable;
import java.net.URL;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * <code>ImageMap</code> represents a collection of image map
 * hotspots. It does not represent the raster image itself.
 *
 * @author Keith James
 * @author Greg Cox
 * @since 1.3
 */
public interface ImageMap
{
    /**
     * <code>RECT</code> indicates a rectangular image map hotspot.
     */
    public static final String RECT = "rect";

    /**
     * <code>CIRCLE</code> indicates a circular image map hotspot.
     */
    public static final String CIRCLE = "circle";

    /**
     * <code>POLY</code> indicates a polygonal image map hotspot.
     */
    public static final String POLY = "poly";

    /**
     * <code>addHotSpot</code> adds a hotspot to the map.
     *
     * @param hotSpot a <code>HotSpot</code>.
     */
    public void addHotSpot(HotSpot hotSpot);

    /**
     * <code>hotSpots</code> iterates over the hotspots in the map
     *
     * @return an <code>Iterator</code>.
     */
    public Iterator hotSpots();

    /**
     * <p><code>HotSpot</code>s represent an image map hotspot. For
     * example (in server-side map format):</p>
     *
     * <p>rect http://www.biojava.org 0,0 100,20</p>
     *
     * <p>A user object may be set for each hot spot. This would
     * typically contain extra data used to construct a representation
     * of the hotspot in a document or application. For example, in an
     * image map representing Blast search results the user object
     * could be a sequence in a database. In an HTML document the user
     * object could be used to assign values to actions such as
     * mouseover.</p>
     */
    public static final class HotSpot implements Serializable
    {
        private String type;
        private URL url;
        private Integer [] coordinates;
        private Object userObject;

        /**
         * Creates a new <code>HotSpot</code> with a null user object.
         *
         * @param type a <code>String</code> of hotspot. The only
         * valid arguments are ImageMap.RECT, ImageMap.CIRCLE or
         * ImageMap.POLY (checked by object reference equalty);
         * @param url a <code>URL</code> target.
         * @param coordinates an <code>Integer []</code> array of
         * hotspot coordinates, in order.
         */
        public HotSpot(String type, URL url, Integer [] coordinates)
        {
            if (! (type == RECT || type == CIRCLE || type == POLY))
                throw new IllegalArgumentException("Failed to create HotSpot. Constructor was passed an invalid type '"
                                                   + type + "'");

            if (! (coordinates.length % 2 == 0))
                throw new IllegalArgumentException("Failed to create HotSpot. The coordinates array contained an odd number of points");

            this.type = type;
            this.url = url;
            this.coordinates = coordinates;
        }

        /**
         * Creates a new <code>HotSpot</code>.
         *
         * @param type a <code>String</code> of hotspot. The only
         * valid arguments are ImageMap.RECT, ImageMap.CIRCLE or
         * ImageMap.POLY (checked by object reference equalty);
         * @param url a <code>URL</code> target.
         * @param coordinates an <code>Integer []</code> array of
         * hotspot coordinates, in order.
         * @param userObject an <code>Object</code>
         */
        public HotSpot(String type, URL url, Integer [] coordinates, Object userObject)
        {
            this(type, url, coordinates);
            this.userObject = userObject;
        }

        /**
         * <code>getType</code> returns the type of hotspot.
         *
         * @return a <code>String</code>.
         */
        public String getType()
        {
            return type;
        }

        /**
         * <code>getURL</code> returns the hotspot URL.
         *
         * @return a <code>URL</code>.
         */
        public URL getURL()
        {
            return url;
        }

        /**
         * <code>getCoordinates</code> returns the hotspot coordinates.
         *
         * @return an <code>Integer []</code> array.
         */
        public Integer [] getCoordinates()
        {
            return coordinates;
        }

        /**
         * <code>getUserObject</code> returns the current user object
         * (or null).
         *
         * @return an <code>Object</code>.
         */
        public Object getUserObject()
        {
            return userObject;
        }

        /**
         * <code>setUserObject</code> sets the user object.
         *
         * @param userObject an <code>Object</code>.
         */
        public void setUserObject(Object userObject)
        {
            this.userObject = userObject;
        }

        public String toString()
        {
            return  "HotSpot to " + url.toString();
        }
    }

    /**
     * <code>ClientSide</code> represents a client-side style image
     * map.
     */
    public static class ClientSide implements ImageMap, Serializable
    {
        private String name;
        private List hotSpots;

        /**
         * Creates a new <code>ClientSide</code> image map.
         *
         * @param name a <code>String</code> name by which the map
         * will be known.
         */
        public ClientSide(String name)
        {
            this.name = name;
            hotSpots = new ArrayList();
        }

        public void addHotSpot(HotSpot hotSpot)
        {
            hotSpots.add(hotSpot);
        }

        public Iterator hotSpots()
        {
            return hotSpots.iterator();
        }

        public String toString()
        {
            StringBuffer sb = new StringBuffer();

            sb.append("<map name=");
            sb.append("\"");
            sb.append(name);
            sb.append("\">\n");

            for (Iterator hi = hotSpots.iterator(); hi.hasNext();)
            {
                HotSpot hs = (HotSpot) hi.next();
                Integer [] coords = hs.getCoordinates();

                sb.append("<area shape=\"");
                sb.append(hs.getType());
                sb.append("\" href=\"");
                sb.append(hs.getURL().toString());
                sb.append("\" coords=\"");

                int lastDelim = coords.length - 1;
                char delim = ',';

                for (int i = 0; i < coords.length; i++)
                {
                    sb.append(coords[i]);

                    if (i <= lastDelim)
                        sb.append(delim);
                }
                sb.append("\">\n");
            }

            sb.append("</map>");

            return sb.substring(0);
        }
    }

    /**
     * <code>ServerSide</code> represents a server-side style image
     * map.
     */
    public static class ServerSide implements ImageMap, Serializable
    {
        private List hotSpots;

        /**
         * Creates a new <code>ServerSide</code> image map.
         */
        public ServerSide()
        {
            hotSpots = new ArrayList();
        }

        public void addHotSpot(HotSpot hotSpot)
        {
            hotSpots.add(hotSpot);
        }

        public Iterator hotSpots()
        {
            return hotSpots.iterator();
        }

        public String toString()
        {
            StringBuffer sb = new StringBuffer();

            for (Iterator hi = hotSpots.iterator(); hi.hasNext();)
            {
                HotSpot hs = (HotSpot) hi.next();
                Integer [] coords = hs.getCoordinates();

                sb.append(hs.getType());
                sb.append(" ");
                sb.append(hs.getURL().toString());
                sb.append(" ");

                for (int i = 0; i < (coords.length - 1); i += 2)
                {
                    sb.append(coords[i]);
                    sb.append(",");
                    sb.append(coords[i + 1]);
                    sb.append(" ");
                }
            }

            return sb.substring(0).trim();
        }
    }
}
