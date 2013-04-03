// $Id: Distribution.java,v 1.17 2009/10/26 23:29:39 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// Copyright (C) 2000-2001 Washington University School of Medicine
// and Howard Hughes Medical Institute
// All rights reserved
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: cmzmasek@yahoo.com
// WWW: www.phylosoft.org/forester

package org.forester.phylogeny.data;

import java.io.IOException;
import java.io.Writer;
import java.math.BigDecimal;

import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.util.ForesterUtil;

public class Distribution implements PhylogenyData {

    private String     _desc;
    private BigDecimal _latitude;
    private BigDecimal _longitude;
    private BigDecimal _altitude;
    private String     _geodetic_datum;

    public Distribution( final String desc ) {
        _desc = desc;
        _latitude = null;
        _longitude = null;
        _altitude = null;
        _geodetic_datum = "";
    }

    public Distribution( final String desc,
                         final BigDecimal latitude,
                         final BigDecimal longitude,
                         final BigDecimal altitude,
                         final String geodetic_datum ) {
        _desc = desc;
        _latitude = latitude;
        _longitude = longitude;
        _altitude = altitude;
        _geodetic_datum = geodetic_datum;
    }

    @Override
    public StringBuffer asSimpleText() {
        if ( ( getLatitude() != null ) || ( getLongitude() != null ) || ( getAltitude() != null ) ) {
            return new StringBuffer( getDesc() + " [" + getLatitude().toPlainString() + ","
                    + getLongitude().toPlainString() + "," + getAltitude().toPlainString() + "] [" + getGeodeticDatum()
                    + "]" );
        }
        else {
            return new StringBuffer( getDesc() );
        }
    }

    @Override
    public StringBuffer asText() {
        return asSimpleText();
    }

    @Override
    public PhylogenyData copy() {
        return new Distribution( new String( getDesc() ),
                                 new BigDecimal( getLatitude().toPlainString() ),
                                 new BigDecimal( getLongitude().toPlainString() ),
                                 new BigDecimal( getAltitude().toPlainString() ),
                                 new String( getGeodeticDatum() ) );
    }

    public BigDecimal getAltitude() {
        return _altitude;
    }

    public String getDesc() {
        return _desc;
    }

    public String getGeodeticDatum() {
        return _geodetic_datum;
    }

    public BigDecimal getLatitude() {
        return _latitude;
    }

    public BigDecimal getLongitude() {
        return _longitude;
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        throw new UnsupportedOperationException();
    }

    public void setAltitude( final BigDecimal altitude ) {
        _altitude = altitude;
    }

    public void setDescription( final String desc ) {
        _desc = desc;
    }

    public void setGeodeticDatum( final String datum ) {
        _geodetic_datum = datum;
    }

    public void setLatitude( final BigDecimal latitude ) {
        _latitude = latitude;
    }

    public void setLongitude( final BigDecimal longitud ) {
        _longitude = longitud;
    }

    @Override
    public StringBuffer toNHX() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendOpen( writer, PhyloXmlMapping.DISTRIBUTION );
        if ( !ForesterUtil.isEmpty( getDesc() ) ) {
            PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.DISTRIBUTION_DESC, getDesc(), indentation );
        }
        if ( ( getLatitude() != null ) || ( getLongitude() != null ) || ( getAltitude() != null ) ) {
            PhylogenyDataUtil.appendOpen( writer,
                                          PhyloXmlMapping.POINT,
                                          PhyloXmlMapping.POINT_GEODETIC_DATUM,
                                          getGeodeticDatum() );
            PhylogenyDataUtil.appendElement( writer,
                                             PhyloXmlMapping.POINT_LATITUDE,
                                             getLatitude().toPlainString(),
                                             indentation );
            PhylogenyDataUtil.appendElement( writer,
                                             PhyloXmlMapping.POINT_LONGITUDE,
                                             getLongitude().toPlainString(),
                                             indentation );
            PhylogenyDataUtil.appendElement( writer,
                                             PhyloXmlMapping.POINT_ALTITUDE,
                                             getAltitude().toPlainString(),
                                             indentation );
            PhylogenyDataUtil.appendClose( writer, PhyloXmlMapping.POINT );
        }
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendClose( writer, PhyloXmlMapping.DISTRIBUTION );
    }

    @Override
    public String toString() {
        return asSimpleText().toString();
    }
}
