/* -*- c++ -*- */
/*
 * Copyright 2020 Free Software Foundation, Inc.
 *
 * This file is part of GNU Radio
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tcp_sink_impl.h"
#include <gnuradio/io_signature.h>

#include <boost/format.hpp>
#include <chrono>
#include <sstream>
#include <thread>

namespace gr {
namespace network {

tcp_sink::sptr tcp_sink::make(
    size_t itemsize, size_t veclen, const std::string& host, int port, int sinkmode)
{
    return gnuradio::make_block_sptr<tcp_sink_impl>(
        itemsize, veclen, host, port, sinkmode);
}

/*
 * The private constructor
 */
tcp_sink_impl::tcp_sink_impl(
    size_t itemsize, size_t veclen, const std::string& host, int port, int sinkmode)
    : gr::sync_block("tcp_sink",
                     gr::io_signature::make(1, 1, itemsize * veclen),
                     gr::io_signature::make(0, 0, 0)),
      d_itemsize(itemsize),
      d_veclen(veclen),
      d_host(host),
      d_port(port),
      d_sinkmode(sinkmode),
      d_thread_running(false),
      d_stop_thread(false),
      d_listener_thread(NULL),
      d_start_new_listener(false),
      d_initial_connection(true)
{
    d_block_size = d_itemsize * d_veclen;

        start_flag1 = 0; s_cycle1=0;
    freq_flag1 = 0;
    meanVec_flag1 = 0;
    ch_flag1 = 0;
    oneCount1 = 0;
    chCounter1 = -1;
    chCounter2 = 0;
    meanVecCount1=0;

    local_meanVec1 = (gr_complex*) malloc (sizeof(gr_complex)*16384);
    local_ch_idx1 = (gr_complex*) malloc (sizeof(gr_complex)*16384);
    local_ch_cen1 = (gr_complex*) malloc (sizeof(gr_complex)*16384);
    local_ch_mag1 = (gr_complex*) malloc (sizeof(gr_complex)*16384);
    local_ch_bw1 = (gr_complex*) malloc (sizeof(gr_complex)*16384);

    local_meanVec01 = (float*) malloc (sizeof(float)*16384);
    for(int i=0;i<=16383;++i) local_meanVec01[i]=0;

    sumVec = (float*) malloc (sizeof(float)*128);
    for(int i=0;i<=128-1;++i) sumVec[i]=0;

    masBuffer = (float*) malloc (sizeof(float)*212*32);
    for(int i=0;i<=212*32-1;++i) masBuffer[i]=0;

    freqCount=0;
    tempCount=wfCount=tempSum=0;
    buffCount=0;
    firstTimeFlag=0;
    maxBuffCount=20000;
    readCount=-1;

}

bool tcp_sink_impl::start()
{
    if (d_sinkmode == TCPSINKMODE_CLIENT) {
        // In this mode, we're connecting to a remote TCP service listener
        // as a client.
        std::stringstream msg;

        msg << "[TCP Sink] connecting to " << d_host << " on port " << d_port;
        GR_LOG_INFO(d_logger, msg.str());

        boost::system::error_code err;
        d_tcpsocket = new boost::asio::ip::tcp::socket(d_io_service);

        std::string s_port = (boost::format("%d") % d_port).str();
        boost::asio::ip::tcp::resolver resolver(d_io_service);
        boost::asio::ip::tcp::resolver::query query(
            d_host, s_port, boost::asio::ip::resolver_query_base::passive);

        d_endpoint = *resolver.resolve(query, err);

        if (err) {
            throw std::runtime_error(
                std::string("[TCP Sink] Unable to resolve host/IP: ") + err.message());
        }

        if (d_host.find(":") != std::string::npos)
            d_is_ipv6 = true;
        else {
            // This block supports a check that a name rather than an IP is provided.
            // the endpoint is then checked after the resolver is done.
            if (d_endpoint.address().is_v6())
                d_is_ipv6 = true;
            else
                d_is_ipv6 = false;
        }

        d_tcpsocket->connect(d_endpoint, err);
        if (err) {
            throw std::runtime_error(std::string("[TCP Sink] Connection error: ") +
                                     err.message());
        }

        d_connected = true;

        boost::asio::socket_base::keep_alive option(true);
        d_tcpsocket->set_option(option);
    } else {
        // In this mode, we're starting a local port listener and waiting
        // for inbound connections.
        d_start_new_listener = true;
        d_listener_thread = new boost::thread([this] { run_listener(); });
    }

    return true;
}

void tcp_sink_impl::run_listener()
{
    d_thread_running = true;

    while (!d_stop_thread) {
        // this will block
        if (d_start_new_listener) {
            d_start_new_listener = false;
            connect(d_initial_connection);
            d_initial_connection = false;
        } else
            std::this_thread::sleep_for(std::chrono::microseconds(10));
    }

    d_thread_running = false;
}

void tcp_sink_impl::accept_handler(boost::asio::ip::tcp::socket* new_connection,
                                   const boost::system::error_code& error)
{
    if (!error) {
        GR_LOG_INFO(d_logger, "Client connection received.");

        // Accept succeeded.
        d_tcpsocket = new_connection;

        boost::asio::socket_base::keep_alive option(true);
        d_tcpsocket->set_option(option);
        d_connected = true;

    } else {
        std::stringstream msg;
        msg << "Error code " << error << " accepting TCP session.";
        GR_LOG_ERROR(d_logger, msg.str());

        // Boost made a copy so we have to clean up
        delete new_connection;

        // safety settings.
        d_connected = false;
        d_tcpsocket = NULL;
    }
}

void tcp_sink_impl::connect(bool initial_connection)
{
    std::stringstream msg;
    msg << "Waiting for connection on port " << d_port;
    GR_LOG_INFO(d_logger, msg.str());

    if (initial_connection) {
        if (d_is_ipv6)
            d_acceptor = new boost::asio::ip::tcp::acceptor(
                d_io_service,
                boost::asio::ip::tcp::endpoint(boost::asio::ip::tcp::v6(), d_port));
        else
            d_acceptor = new boost::asio::ip::tcp::acceptor(
                d_io_service,
                boost::asio::ip::tcp::endpoint(boost::asio::ip::tcp::v4(), d_port));
    } else {
        d_io_service.reset();
    }

    if (d_tcpsocket) {
        delete d_tcpsocket;
    }
    d_tcpsocket = NULL;
    d_connected = false;

    boost::asio::ip::tcp::socket* tmpSocket =
        new boost::asio::ip::tcp::socket(d_io_service);
    d_acceptor->async_accept(*tmpSocket,
                             boost::bind(&tcp_sink_impl::accept_handler,
                                         this,
                                         tmpSocket,
                                         boost::asio::placeholders::error));

    d_io_service.run();
}

/*
 * Our virtual destructor.
 */
tcp_sink_impl::~tcp_sink_impl() { stop();

    free (local_meanVec1);
    free (local_ch_idx1);
    free (local_ch_cen1);
    free (local_ch_mag1);
    free (local_ch_bw1);
    free (local_meanVec01);
    free (sumVec);
    free (masBuffer);

}

bool tcp_sink_impl::stop()
{
    if (d_thread_running) {
        d_stop_thread = true;
    }

    if (d_tcpsocket) {
        d_tcpsocket->close();
        delete d_tcpsocket;
        d_tcpsocket = NULL;
    }

    d_io_service.reset();
    d_io_service.stop();

    if (d_acceptor) {
        delete d_acceptor;
        d_acceptor = NULL;
    }

    if (d_listener_thread) {
        while (d_thread_running)
            std::this_thread::sleep_for(std::chrono::microseconds(5));

        delete d_listener_thread;
        d_listener_thread = NULL;
    }

    return true;
}

int tcp_sink_impl::work(int noutput_items,
                        gr_vector_const_void_star& input_items,
                        gr_vector_void_star& output_items)
{
    gr::thread::scoped_lock guard(d_setlock);



    const gr_complex* in = (const gr_complex*)input_items[0];

    for (int i = 0; i <= noutput_items-1; ++i){

           if(start_flag1==0 && in[i]==gr_complex(1,0)) {
        	start_flag1=1;
        	oneCount1++;    //  std::cout<<s_cycle1+i<<" first "<<i<<"  "<<noutput_items<<"  "<<in[i-1]<<"  "<<in[i]<<"  "<<in[i+1]<<"\n";
        }

        else if(start_flag1==1 && in[i]==gr_complex(1,0) && ch_flag1 == 0) {
        	oneCount1++;
        	if(oneCount1==16383){
						//	 std::cout<<s_cycle1+i<<" second "<<i<<"  "<<noutput_items<<"  "<<in[i-1]<<"  "<<in[i]<<"  "<<in[i+1]<<"\n";
                freq_flag1=1;
        	}
        }

        else if(freq_flag1==1 && start_flag1==1){
            local_freq1 = in[i];  	//  std::cout<<s_cycle1+i<<" third "<<i<<"  "<<noutput_items<<"  "<<local_freq1<<"  "<<in[i-1]<<"  "<<in[i]<<"  "<<in[i+1]<<"\n";
        	meanVec_flag1=1;
        	freq_flag1=0;
        	if((int)local_freq1.real() == 40000000) {
                freqCount=0;

                buffCount=0;

            }
        	freqCount++;	// std::cout<<local_freq1<<"   "<<tempCount<<"\n";
       //     buffCount = (freqCount-1)*212;
            masBuffer[buffCount]=111111; buffCount++;
            masBuffer[buffCount]=local_freq1.real(); buffCount++;
        }

        else if(meanVec_flag1==1 && start_flag1==1){
        	local_meanVec1[meanVecCount1]=in[i];
        	local_meanVec01[meanVecCount1]=in[i].real();
        	meanVecCount1++;
        	tempSum= tempSum+local_meanVec01[meanVecCount1]; tempCount++;
        	if(tempCount==128){
        		tempSum/=128.0;
        		sumVec[wfCount] = 100*tempSum;	//  std::cout<<meanVecCount1<<"   "<<wfCount<<"\n";
                masBuffer[buffCount] = sumVec[wfCount]; buffCount++;
        		tempCount=tempSum=0;
        		wfCount++;
        		if(wfCount==128-1) { wfCount=0; }
        	}

        	if(meanVecCount1==16384){	//  std::cout<<s_cycle1+i<<" fourth "<<i<<"  "<<noutput_items<<"  "<<in[i-1]<<"  "<<in[i]<<"  "<<in[i+1]<<"\n";
        		ch_flag1=1;
        		meanVec_flag1=0;
        	}
        }

        else if (ch_flag1 == 1 && start_flag1 == 1 && chCounter1 == -1){
            total_chs1 = in[i].real();		//  std::cout<<s_cycle1+i<<" fifth "<<i<<"  "<<total_chs1<<"  "<<in[i-1]<<"  "<<in[i]<<"  "<<in[i+1]<<"\n";
            masBuffer[buffCount] = total_chs1; buffCount++;
            chCounter1++;
            if(total_chs1==0){			//  std::cout<<total_chs1<<" dummy_fifth "<<local_freq1<<"\n";
                ch_flag1 = 0;
                start_flag1 = 0;
                oneCount1 = 0;
                meanVecCount1= 0;
                chCounter1 = -1;
                local_freq1 = -1;

                if(freqCount==32){
                    masBuffer[buffCount]=999999; buffCount++;
                    maxBuffCount=buffCount;
                }
                else { masBuffer[buffCount]=222222; buffCount++; }

            }
        }

        else if (ch_flag1 == 1 && start_flag1 == 1 && chCounter1 >= 0){
            if(chCounter2==0) { local_ch_idx1[chCounter1] = in[i]; chCounter2++;}
            else if(chCounter2==1) { local_ch_cen1[chCounter1] = in[i]; chCounter2++;}
            else if(chCounter2==2) { local_ch_mag1[chCounter1] = in[i]; chCounter2++;}
            else if(chCounter2==3) {
                local_ch_bw1[chCounter1] = in[i];

         //       std::cout<<"serial:"<<total_chs1<<" dummy_sixth "<<local_freq1<<"  "<<"idx:"<<local_ch_idx1[chCounter1]<<" center freq: "<<local_ch_cen1[chCounter1]<<" mag: "<<local_ch_mag1[chCounter1]<<" bw: "<<local_ch_bw1[chCounter1] <<"\n";
                if(chCounter1<=19){
                    masBuffer[buffCount]=local_ch_idx1[chCounter1].real(); buffCount++;
                    masBuffer[buffCount]=local_ch_cen1[chCounter1].real(); buffCount++;
                    masBuffer[buffCount]=local_ch_mag1[chCounter1].real(); buffCount++;
                    masBuffer[buffCount]=local_ch_bw1[chCounter1].real(); buffCount++;
                }

                chCounter2 = 0;
                chCounter1++;
                if(chCounter1 == total_chs1){		// std::cout<<s_cycle1+i<<" sixth "<<i<<"  "<<noutput_items<<"  "<<chCounter1<<"  "<<in[i-1]<<"  "<<in[i]<<"  "<<in[i+1]<<"\n";
                    ch_flag1 = 0;
                    start_flag1 = 0;
                    oneCount1 = 0;
                    meanVecCount1= 0;
                    chCounter1 = -1;
                    local_freq1 = -1;

                    if(freqCount==32) {
                        masBuffer[buffCount]=999999; buffCount++;
                        maxBuffCount=buffCount;
                        // std::cout<<buffCount<<"  "<<masBuffer[buffCount-1]<<"  "<<masBuffer[buffCount-100]<<"\n";
                       // char* p1; p1 = (char*) &masBuffer[buffCount]; printf("%x",*(p1+2)); std::cout<<"\n";
                    }
                    else { masBuffer[buffCount]=222222; buffCount++; }
                }
            }

        }

    }
     s_cycle1+= noutput_items;



    if (!d_connected)
        return noutput_items;

    if(buffCount>5 && buffCount<213 && firstTimeFlag==0){
        readCount=0;
        maxBuffCount=20000;
        firstTimeFlag=1;
      }

    int bytes_written;
    ec.clear();

    char* p_buff; p_buff = (char*) &masBuffer[readCount];

   // std::cout<<readCount<<"  "<<buffCount<<"\n";


    if (!ec && readCount>=0) {


            int delta, topLimit;
            if(maxBuffCount==20000) topLimit=buffCount;
            else topLimit=maxBuffCount;

            delta=topLimit-readCount;
            if(delta>5) delta=5;
      //      else if(topLimit==buffCount) delta=delta+1;

            if(delta>0){
                bytes_written = boost::asio::write(*d_tcpsocket, boost::asio::buffer((const void*)p_buff, delta*4), ec);
            //    std::cout<<delta<<"  "<<masBuffer[readCount]<<"  "<<masBuffer[readCount+1]<<"  "<<masBuffer[readCount+2]<<"  "<<masBuffer[readCount+3]<<"  "<<masBuffer[readCount+4]<<"  "<<masBuffer[readCount+5]<<" \n ";
                readCount=readCount+delta;
            }
            if((readCount>=maxBuffCount && buffCount==maxBuffCount)|| readCount>=maxBuffCount) {
                firstTimeFlag=0;
                readCount=-1;
            }


      //  std::cout<<"inside: "<<readCount<<"  "<<buffCount<< " "<<maxBuffCount<<"\n";

        if (ec == boost::asio::error::connection_reset ||
            ec == boost::asio::error::broken_pipe) {

            // Connection was reset
                d_connected = false;
                GR_LOG_INFO(d_logger, "Client disconnected. Waiting for new connection.");
                // start waiting for another connection
                d_start_new_listener = true;
        }
    }

    return noutput_items;
}
} /* namespace network */
} /* namespace gr */
