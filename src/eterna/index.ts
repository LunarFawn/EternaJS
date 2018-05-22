import * as log from "loglevel";
import {EternaApp} from "./EternaApp";

const isProduction = process.env.NODE_ENV === 'production';
log.setLevel(isProduction ? 'warn' : 'trace');

let app: EternaApp = new EternaApp();
app.run();