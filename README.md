# Vector Quantization For Neural Signals

Short Description: </br>
In the thesis, we explore the potential of vector quantization (VQ) combined with discrete cosine transform (DCT) and Huffman coding for data reduction in a wireless brain monitoring devices. We look at extracellular tetrode recodrings containing Local Field Potentials (LFP), and Extracellular Action Potentials (EAP), and focus on the compression of the EAP signal as it is the most challenging to transmit/store due to its bandwidth and resolution requirements. We investigate how the proposed compression algorithm affects spike sorting at different compression rates. This GitHub repository contains all of the code needed to reproduce the experiments that have been conducted. 

Dataset used: </br>
The HC-1 dataset has been used for all testing, and can be found at: https://crcns.org/data-sets/hc/hc-1


Thesis Abstract: </br>
Brain disorders like Dementia and Parkinson's affect millions and constitute one of the biggest health problems that the world faces today. Yet, we are limited by a lack of mechanical understanding of the brain . To address this, research is being conducted to develop wireless brain monitoring devices for use on rats. These devices require sending high amounts of data while maintaining a small design, resulting in challenges related to data transmission and power consumption. This thesis explores the compression of extracellular action potentials (EAPs) as a means to mitigate these issues. It presents a novel method that can learn the underlying structure of an EAP signal, so as to compress spikes with higher resolution than the background signal without the need for local spike detection or curation. 

The method leverages the discrete cosine transform (DCT), vector quantization (VQ), and Huffman coding to compress the extracellular tetrode signal of the publicly available in-vivo HC-1 (d533101) dataset. The results show an increase in spike sorting recall from 92.7% to 94% compared to the uncompressed signal while achieving a compression rate (CR) of 32. For higher CRs of 73 and 131, the method maintains spike sorting accuracies of 98% and 83%, respectively, compared to the uncompressed signal. This indicates the method's superiority over existing algorithms that do not use local spike detection. As a next step the thesis recommends testing different transforms, that make the spikes and background signal more separable, to enhance spike sorting at high compression rates and potentially outperform current compression methods that rely on local spike detection. The findings presented in this thesis indicate the potential of such an approach in the field of EAP compression. 

