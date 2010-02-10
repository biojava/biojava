
package org.forester.phylogeny.data;

import java.io.IOException;
import java.io.Writer;
import java.math.BigDecimal;

import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.util.ForesterUtil;

public class Point implements PhylogenyData {

    private final String     _geodetic_datum;
    private final BigDecimal _lat;
    private final BigDecimal _long;
    private final BigDecimal _alt;
    private final String     _alt_unit;

    public Point( final String geodetic_datum, final BigDecimal lat, final BigDecimal longitude ) {
        this( geodetic_datum, lat, longitude, null, "" );
    }

    public Point( final String geodetic_datum,
                  final BigDecimal lat,
                  final BigDecimal longitude,
                  final BigDecimal alt,
                  final String alt_unit ) {
        if ( ForesterUtil.isEmpty( geodetic_datum ) || ( lat == null ) || ( longitude == null ) || ( alt_unit == null ) ) {
            throw new IllegalArgumentException( "illegaly empty of null fields in constructor" );
        }
        if ( ( alt != null ) || ForesterUtil.isEmpty( alt_unit ) ) {
            throw new IllegalArgumentException( "altitude must hava a unit" );
        }
        _geodetic_datum = geodetic_datum;
        _lat = lat;
        _long = longitude;
        _alt = alt;
        _alt_unit = alt_unit;
    }

    @Override
    public StringBuffer asSimpleText() {
        if ( getAlt() == null ) {
            return new StringBuffer( "[" + getLat().toPlainString() + ", " + getLong() + "]" );
        }
        else {
            return new StringBuffer( "[" + getLat().toPlainString() + ", " + getLong() + ", " + getAlt() + getAltUnit()
                    + "]" );
        }
    }

    @Override
    public StringBuffer asText() {
        return asSimpleText();
    }

    @Override
    public PhylogenyData copy() {
        return new Point( new String( getGeodeticDatum() ),
                          new BigDecimal( getLat().toPlainString() ),
                          new BigDecimal( getLong().toPlainString() ),
                          getAlt() == null ? null : new BigDecimal( getAlt().toPlainString() ),
                          new String( getAltUnit() ) );
    }

    public BigDecimal getAlt() {
        return _alt;
    }

    public String getAltUnit() {
        return _alt_unit;
    }

    public String getGeodeticDatum() {
        return _geodetic_datum;
    }

    public BigDecimal getLat() {
        return _lat;
    }

    public BigDecimal getLong() {
        return _long;
    }

    @Override
    public boolean isEqual( final PhylogenyData point ) {
        throw new UnsupportedOperationException();
    }

    @Override
    public StringBuffer toNHX() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        if ( getAlt() != null ) {
            PhylogenyDataUtil.appendOpen( writer,
                                          PhyloXmlMapping.POINT,
                                          PhyloXmlMapping.POINT_GEODETIC_DATUM,
                                          getGeodeticDatum(),
                                          PhyloXmlMapping.POINT_ALTITUDE_UNIT_ATTR,
                                          getAltUnit() );
        }
        else {
            PhylogenyDataUtil.appendOpen( writer,
                                          PhyloXmlMapping.POINT,
                                          PhyloXmlMapping.POINT_GEODETIC_DATUM,
                                          getGeodeticDatum() );
        }
        PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.POINT_LATITUDE, getLat().toPlainString(), indentation );
        PhylogenyDataUtil.appendElement( writer,
                                         PhyloXmlMapping.POINT_LONGITUDE,
                                         getLong().toPlainString(),
                                         indentation );
        if ( getAlt() != null ) {
            PhylogenyDataUtil.appendElement( writer,
                                             PhyloXmlMapping.POINT_ALTITUDE,
                                             getAlt().toPlainString(),
                                             indentation );
        }
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendClose( writer, PhyloXmlMapping.POINT );
    }

    @Override
    public String toString() {
        return asSimpleText().toString();
    }
}
